import io
from contextlib import redirect_stdout
from typing import Optional
from IPython.core.interactiveshell import InteractiveShell #just for the type annotation
from IPython.core.magic import register_line_magic

@register_line_magic # so that we can optionally use get_magic_out as a magic command
def get_magic_outb(command: str, ipy: Optional[InteractiveShell]=None) -> str:
    """ Redirects the stdout of a ipython magic command to a temporary buffer
    Args:
        command: The command to be executed, as if it were in the shell.
            ie `get_magic_out('history 5')` would capture the output of the command
            ```
            In [1]: %history 5
            ```
         ipy (optional): The ipython shell to run the command in. Defaults to getting
            the current shell (by calling `get_ipython` which is bound to the current shell)
            
     Returns:
        The output of %{command}
    """
    ipy = ipy or get_ipython()
    out = io.BytesIO()
    
    with redirect_stdout(out):
        ipy.magic(command)

    return out.getvalue()
