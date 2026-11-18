"""Launch the desktop application with this checkout's portable catalog."""
import multiprocessing
import sys
from pathlib import Path

if not getattr(sys, "frozen", False):
    repository = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repository / "src"))

if __name__ == "__main__":
    multiprocessing.freeze_support()
    arguments = sys.argv[1:]
    if arguments and arguments[0] == '--raspa-worker':
        from cof_vision.cofdraw_raspa_worker import main as worker_main
        raise SystemExit(worker_main(arguments[1:]))
    from cof_vision.cofdraw_desktop import main
    if not getattr(sys, "frozen", False) and not any(
        arg == "--data-dir" or arg.startswith("--data-dir=") for arg in arguments
    ):
        arguments = ["--data-dir", str(repository / "data"), *arguments]
    raise SystemExit(main(arguments))
