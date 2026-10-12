"""OrthoPPInterface: Master CLI dispatcher to pipeline/run_pipeline.py."""
import os
import sys
import subprocess

if __name__ == "__main__":
    pipeline_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
    script = os.path.join(pipeline_dir, "run_pipeline.py")
    cmd = [sys.executable, script] + sys.argv[1:]
    sys.exit(subprocess.call(cmd, cwd=pipeline_dir))
