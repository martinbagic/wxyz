import argparse
import yaml
import pathlib
import json
import itertools

# get filename from cmd
parser = argparse.ArgumentParser("Summoner")
parser.add_argument("filename", type=str)
filename = parser.parse_args().filename

# read yml
path = pathlib.Path(__file__).absolute().parent / filename
with open(path, "r") as f:
    yml = yaml.safe_load(f)

# initialize container for commands
commands = []
command_pattern = (
    "python3 main.py {paramjson} {experiment_name} {experiment_dir} {c} {r}"
)

# collect commands
for line in yml:
    values = list(line["params"].values())
    keys = list(line["params"].keys())
    param_tupls = itertools.product(*values)

    # make command for every parameter combination
    for experiment_i, tupl in enumerate(param_tupls):
        command = command_pattern.format(
            paramjson="'" + json.dumps(dict(zip(keys, tupl))) + "'",
            experiment_name=line["experiment_name"] + "_" + str(experiment_i),
            experiment_dir=line["experiment_dir"],
            c=f"-c {line['c']}" if line["c"] else "",
            r=f"-r {line['r']}" if line["r"] else "",
        )
        commands.append(command)

# write commands to file
with open((path.parent / filename).with_suffix(".summon"), "w") as f:
    f.write("\n".join(commands))

