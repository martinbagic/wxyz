import yaml

def lister(path):
    
    with open(path, "r") as ofile:
        listing = yaml.safe_load(ofile)
        
    print(listing)
    
    for params in listing[0]["extra_params"]:
        print(eval(params[1]))
    
    
lister("xhelp/listing.yml")