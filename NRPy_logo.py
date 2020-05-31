def print_logo(print_to_stdout=True):

    logo_str = """
ooooo      ooo ooooooooo.   ooooooooo.                 88
`888b.     `8' `888   `Y88. `888   `Y88.             888888
 8 `88b.    8   888   .d88'  888   .d88' oooo    ooo   88
 8   `88b.  8   888ooo88P'   888ooo88P'   `88.  .8'
 8     `88b.8   888`88b.     888           `88..8'
 8       `888   888  `88b.   888            `888'
o8o        `8  o888o  o888o o888o            .8'
                                         .o..P'
  NRPy+: Python-based Code Generation    `Y8P'
   for Numerical Relativity... and Beyond!
 - homepage: http://blackholesathome.net
 - download: https://github.com/zachetienne/nrpytutorial
"""
    if print_to_stdout==True:
        print(logo_str)
    else:
        return logo_str
