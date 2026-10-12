"""OrthoPPInterface: Master CLI dispatcher to pipeline/run_pipeline.py."""
import os
import sys
import subprocess

if __name__ == "__main__":
    root_dir = os.path.dirname(os.path.abspath(__file__))
    pipeline_dir = os.path.join(root_dir, "pipeline")
    script = os.path.join(pipeline_dir, "run_pipeline.py")

    args = list(sys.argv[1:])
    for i, arg in enumerate(args):
        if arg == "--config" and i + 1 < len(args):
            val = args[i + 1]
            if os.path.exists(val):
                args[i + 1] = os.path.abspath(val)
            elif os.path.exists(os.path.join(root_dir, val)):
                args[i + 1] = os.path.abspath(os.path.join(root_dir, val))
            elif os.path.exists(os.path.join(pipeline_dir, val)):
                args[i + 1] = os.path.abspath(os.path.join(pipeline_dir, val))
        elif arg.startswith("--config="):
            val = arg.split("=", 1)[1]
            if os.path.exists(val):
                args[i] = f"--config={os.path.abspath(val)}"
            elif os.path.exists(os.path.join(root_dir, val)):
                args[i] = f"--config={os.path.abspath(os.path.join(root_dir, val))}"
            elif os.path.exists(os.path.join(pipeline_dir, val)):
                args[i] = f"--config={os.path.abspath(os.path.join(pipeline_dir, val))}"

    cmd = [sys.executable, script] + args
    sys.exit(subprocess.call(cmd, cwd=pipeline_dir))

