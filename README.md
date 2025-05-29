Repository containing the coordinate data and notebook for robot-assisted stereoelectroencephalography (SEEG) implantations from 2017-2024, totalling 246 subjects and 2658 trajectories. This also includes a subset of 28 patients (289) between 2013-2017 who were implanted in a traditional frame based method, without assistance of the robot.

The project is managed using  [uv](https://github.com/astral-sh/uv), which provides fast dependency resolution and environment management.

To work on the repository, first install `uv`. On macOS and Linux, run:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On Windows, run:
```
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Alternatively, you can install it from PyPI using:

```
pip install uv
```

or with pipx:

```
pipx install uv
```

Once uv is installed, you can clone the repository 

```sh
git clone https://github.com/arunthurai/seeg_accuracy_2025.git
```

Then create and activate a virtual environment:

```
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Can install the package in editable mode:

```
uv pip install -e .
```

Or can run the notebook with this command, once your local copy of this repository.

```
uv run --with jupyter jupyter lab
```