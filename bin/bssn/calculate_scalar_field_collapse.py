from nrpytutorial.dirtree import Dirtree

HOME = 'build/bssn/scalar_field_collapse'


def setup_dirtree(home) -> Dirtree:
    dirtree = Dirtree(home)
    dirtree.add_new_dir('codes', 'root')
    dirtree.add_new_dir('build', 'root')
    dirtree.add_new_dir('bin', 'root')
    dirtree.add_new_dir('output', 'root')
    return dirtree


if __name__ == '__main__':
    dirtree = setup_dirtree(HOME)
    dirtree.create()
