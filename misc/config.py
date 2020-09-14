import yaml

CONFIG = None


class Config:
    def __init__(self, path_default, params_extra):
        with open(path_default, "r") as ifile:
            params_default = yaml.safe_load(ifile)
        params = {**params_default, **params_extra}

        for pkey, pval in params.items():
            setattr(self, pkey, pval)

        assert (
            params["bits_per_locus"] % 2 == 0
        )  # check that you have 2 conceptual chromosomes

        global CONFIG
        CONFIG = self


# def init_params(path_default, params_extra):
#     params = Params(path_default, params_extra)
#     global CONFIG
#     CONFIG = params


# def init(config={}):

#     with open("config.yml", "r") as f:
#         d = yaml.load(f, Loader=yaml.FullLoader)

#     for k, v in config.items():
#         d[k] = v

#     assert d["bits_per_locus"] % 2 == 0  # check that you have 2 conceptual chromosomes

#     return collections.namedtuple("Config", d.keys())(*d.values())


# def init_config(config={}):

#     with open(path / "config.yml", "r") as f:
#         d = yaml.load(f, Loader=yaml.FullLoader)

#     for k, v in config.items():
#         d[k] = v

#     return collections.namedtuple("Config", d.keys())(*d.values())
