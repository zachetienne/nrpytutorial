#------------------------------------------------------------------------------
# Custom whitespace stripping -- Append to bottom of ~/.jupyter/jupyter_notebook_config.py
#  Found here: https://github.com/jupyter/notebook/issues/1455#issuecomment-519891159
#------------------------------------------------------------------------------
def strip_white_space(text):
    return '\n'.join([line.rstrip() for line in text.split('\n')])

def scrub_output_pre_save(model=None, **kwargs):
    """Auto strip trailing white space before saving."""
    # If we are dealing with a notebook #
    if model['type'] == 'notebook':
        # Don't run on anything other than nbformat version 4 #
        if model['content']['nbformat'] != 4:
            print("Skipping white space stripping since `nbformat` != 4.")
            return
        # Apply function to every cell #
        print("Stripping white space on a notebook.")
        for cell in model['content']['cells']:
            if cell['cell_type'] != 'code': continue
            cell['source'] = strip_white_space(cell['source'])
    # If we are dealing with a file #
    if model['type'] == 'file':
        if model['format'] == 'text':
            print("Stripping white space on a file.")
            model['content'] = strip_white_space(model['content'])

c.ContentsManager.pre_save_hook = scrub_output_pre_save
