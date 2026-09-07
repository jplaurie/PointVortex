"""Exercise real solver output/restart behavior; accepts an executable or MPI command."""
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def check_times(actual, expected):
    assert len(actual) == len(expected), (actual, expected)
    assert all(math.isclose(a, b, rel_tol=1e-12, abs_tol=0.0) for a, b in zip(actual, expected)), (actual, expected)


def csv_times(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    seen = set()
    for row in rows:
        assert all(math.isfinite(float(value)) for value in row.values()), row
        identity = (row['time'], row.get('index'))
        assert identity not in seen, ('duplicate output row', row)
        seen.add(identity)
    times = list(dict.fromkeys(float(row['time']) for row in rows))
    assert times == sorted(times)
    return times


def csv_frames(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return list(dict.fromkeys(int(row['frame']) for row in rows))


def checkpoint_times(directory):
    return [float(path.read_text().splitlines()[1].split()[1])
            for path in sorted(directory.glob('checkpoint_*.dat'))]


def main(command):
    with tempfile.TemporaryDirectory(prefix='point_vortex_output_test_') as temporary:
        root = Path(temporary)

        def run(name, extra='', integrator='rk4', success=True):
            directory = root / name
            params = root / (name.replace('/', '_') + '.params')
            params.write_text(
                f'N 2\nboundaryCondition infinite\nintegrator {integrator}\n'
                'timeStep 0.007\nnumThreads 1\nendTime 0.13\noutputTime 0.04\n'
                f'runDirectory {directory}\n' + extra)
            result = subprocess.run([*command, str(params)], text=True, capture_output=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return directory, result

        for integrator in ('rk4', 'dopri5'):
            # No parent directories exist: the solver creates all of them.
            common, result = run(integrator + '_common', integrator=integrator)
            expected = [0, 0.04, 0.08, 0.12, 0.13]
            check_times(csv_times(common / 'trajectory.csv'), expected)
            check_times(csv_times(common / 'diagnostics.csv'), expected)
            check_times(checkpoint_times(common / 'checkpoints'), expected)
            assert str(common / 'trajectory.csv') in result.stdout
            assert str(common / 'diagnostics.csv') in result.stdout
            assert str(common / 'checkpoints') in result.stdout

            split, _ = run(integrator + '_split', 'diagnosticsTime 0.03\ncheckpointTime 0.05\n', integrator)
            check_times(csv_times(split / 'trajectory.csv'), expected)
            check_times(csv_times(split / 'diagnostics.csv'), [0, .03, .06, .09, .12, .13])
            check_times(checkpoint_times(split / 'checkpoints'), [0, .05, .1, .13])
            source = split / 'checkpoints/checkpoint_00000001.dat'
            branch, _ = run(integrator + '_branch',
                            f'diagnosticsTime 0.03\ncheckpointTime 0.05\nrestartFile {source}\n', integrator)
            check_times(csv_times(branch / 'trajectory.csv'), [.05, .08, .12, .13])
            check_times(csv_times(branch / 'diagnostics.csv'), [.05, .06, .09, .12, .13])
            check_times(checkpoint_times(branch / 'checkpoints'), [.1, .13])
            assert (branch / 'checkpoints/checkpoint_00000002.dat').exists()
            with (split / 'trajectory.csv').open() as stream:
                full_rows = list(csv.DictReader(stream))[-2:]
            with (branch / 'trajectory.csv').open() as stream:
                branch_rows = list(csv.DictReader(stream))[-2:]
            for full, resumed in zip(full_rows, branch_rows):
                for key in ('x', 'y', 'circulation', 'u', 'v'):
                    assert math.isclose(float(full[key]), float(resumed[key]), abs_tol=1e-12)

            changed, _ = run(integrator + '_changed',
                             f'outputTime 0.02\ndiagnosticsTime 0.04\ncheckpointTime 0.03\nrestartFile {source}\n', integrator)
            check_times(csv_times(changed / 'trajectory.csv'), [.05, .07, .09, .11, .13])
            check_times(csv_times(changed / 'diagnostics.csv'), [.05, .09, .13])
            check_times(checkpoint_times(changed / 'checkpoints'), [.08, .11, .13])

            # Version 4 has only the shared next-output clock.
            old = root / (integrator + '_v4.dat')
            old.write_text('\n'.join(line.replace('POINT_VORTEX_CHECKPOINT 6', 'POINT_VORTEX_CHECKPOINT 4')
                                     for line in (common / 'checkpoints/checkpoint_00000001.dat').read_text().splitlines()
                                     if not line.startswith(('output_schedule ', 'event_index '))) + '\n')
            legacy, _ = run(integrator + '_legacy', f'restartFile {old}\n', integrator)
            check_times(csv_times(legacy / 'trajectory.csv'), [.04, .08, .12, .13])
            check_times(csv_times(legacy / 'diagnostics.csv'), [.04, .08, .12, .13])
            check_times(checkpoint_times(legacy / 'checkpoints'), [.08, .12, .13])

        for key in ('diagnosticsTime', 'checkpointTime'):
            for value in ('0', '-1', 'nan', 'inf'):
                _, result = run(key + value, f'{key} {value}\n', success=False)
                assert 'error:' in result.stderr

        _, result = run('legacy_output_key', 'outputFile obsolete.csv\n', success=False)
        assert 'unknown parameter: outputFile' in result.stderr

        # Existing managed output must not be replaced unless explicitly requested.
        protected = root / 'protected'
        protected.mkdir()
        sentinel = protected / 'diagnostics.csv'
        sentinel.write_text('previous results\n')
        _, result = run('protected', success=False)
        assert 'already contains solver output' in result.stderr
        assert sentinel.read_text() == 'previous results\n'
        assert not (protected / 'trajectory.csv').exists()
        _, result = run('protected', 'overwriteRun true\n')
        assert result.returncode == 0
        assert sentinel.read_text().startswith('time,frame,')

        malformed = root / 'malformed.dat'
        malformed.write_text('not a vortex\n0 0 1\n')
        _, result = run('malformed', f'initialConditionFile {malformed}\n', success=False)
        assert 'invalid initial condition on line 1' in result.stderr

        # End times differing from a scheduled output only by roundoff need one final frame.
        rounded, _ = run('rounded', 'endTime 0.1200000000000001\n')
        check_times(csv_times(rounded / 'trajectory.csv'), [0, .04, .08, .1200000000000001])
        # Tolerances must scale with time, not with an artificial unit-time floor.
        tiny, _ = run('tiny', 'N 1\ntimeStep 1e-17\nendTime 1e-15\noutputTime 2e-16\n')
        check_times(csv_times(tiny / 'trajectory.csv'), [0, 2e-16, 4e-16, 6e-16, 8e-16, 1e-15])

        for key, value in [('output_index', '-1'), ('event_index', '-1'),
                           ('vortex_count', '999999999999999999'), ('next_output_time', '0'),
                           ('accepted_steps', '-1')]:
            corrupt = root / (key + '.dat')
            lines = source.read_text().splitlines()
            corrupt.write_text('\n'.join(f'{key} {value}' if line.startswith(key + ' ') else line
                                         for line in lines) + '\n')
            run('corrupt_' + key, f'restartFile {corrupt}\n', 'dopri5', success=False)

        initial = root / 'removed.dat'
        initial.write_text('0 0 1\n0.001 0 -1\n')
        removed, _ = run('removed', f'initialConditionFile {initial}\n'
                         'dipoleRemoval true\ndipoleRemovalDistance 0.01\n'
                         'diagnosticsTime 0.03\ncheckpointTime 0.05\n')
        assert csv_times(removed / 'trajectory.csv') == []
        check_times(csv_times(removed / 'diagnostics.csv'), [0, .03, .06, .09, .12, .13])
        empty_source = removed / 'checkpoints/checkpoint_00000001.dat'
        resumed, _ = run('removed_restart', f'restartFile {empty_source}\n'
                         'dipoleRemoval true\ndipoleRemovalDistance 0.01\n'
                         'diagnosticsTime 0.03\ncheckpointTime 0.05\n')
        assert csv_times(resumed / 'trajectory.csv') == []
        check_times(checkpoint_times(resumed / 'checkpoints'), [.1, .13])

        zero, _ = run('zero', 'endTime 0\n')
        check_times(csv_times(zero / 'trajectory.csv'), [0])
        check_times(csv_times(zero / 'diagnostics.csv'), [0])
        check_times(checkpoint_times(zero / 'checkpoints'), [0])

        # Every run owns all outputs and records exactly how it was started.
        managed = root / 'managed_run'
        managed_params = root / 'managed.params'
        managed_params.write_text(
            'N 2\nboundaryCondition infinite\nintegrator rk4\ntimeStep 0.007\n'
            'endTime 0.08\noutputTime 0.04\nnumThreads 1\n'
            f'runDirectory {managed}\n')
        result = subprocess.run([*command, str(managed_params)], text=True, capture_output=True,
                                timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        trajectory = managed / 'trajectory.csv'
        diagnostics = managed / 'diagnostics.csv'
        assert trajectory.exists() and diagnostics.exists()
        assert (managed / 'checkpoints/checkpoint_00000002.dat').exists()
        assert csv_frames(trajectory) == [0, 1, 2], csv_frames(trajectory)
        assert csv_frames(diagnostics) == [0, 1, 2], csv_frames(diagnostics)
        record = (managed / 'resolved_parameters.txt').read_text()
        assert 'POINT_VORTEX_RUN_RECORD 1' in record
        assert 'backend ' in record
        assert 'trajectory_file ' in record and 'checkpoints' in record
        segment = managed / 'segments/segment_00000001/resolved_parameters.txt'
        assert segment.exists()

        managed_branch = root / 'managed_branch'
        branch_params = root / 'managed_branch.params'
        branch_params.write_text(
            'boundaryCondition infinite\nintegrator rk4\ntimeStep 0.007\nendTime 0.12\n'
            'outputTime 0.04\nnumThreads 1\n'
            f'restartFile {managed}/checkpoints/checkpoint_00000001.dat\n'
            f'runDirectory {managed_branch}\n')
        result = subprocess.run([*command, str(branch_params)], text=True, capture_output=True,
                                timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        branch_record = (managed_branch / 'resolved_parameters.txt').read_text()
        assert 'restarting 1' in branch_record
        assert csv_frames(managed_branch / 'trajectory.csv') == [1, 2, 3], \
            csv_frames(managed_branch / 'trajectory.csv')

    print('output scheduling, automatic directories, overwrite protection, and restart tests passed')


if __name__ == '__main__':
    main(sys.argv[1:])
