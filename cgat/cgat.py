'''
cgat.py - Computational Genomics Analysis Tools
===============================================

:Tags: Genomics

To use a specific tool, type::

    cgat <tool> [tool options] [tool arguments]

Tools are grouped by keywords. For this message and a list of
available keywords type::

    cgat --help

For a list of tools matching a certain keyword, type::

   cgat --help <keyword>

or::

   cgat --help all

for a list of all available tools.

To get help for a specific tool, type::

    cgat <tool> --help
'''
import os
import sys
import re
import glob
import importlib.util
import collections
import cgatcore.iotools as iotools
import cgat


def mapKeyword2Script(path):
    '''Collect keywords from scripts.'''

    map_keyword2script = collections.defaultdict(list)

    for script in glob.glob(os.path.join(path, "*.py")):
        s = os.path.basename(script)[:-3]
        with iotools.open_file(script, 'r') as inf:
            data = [x for x in inf.readlines(10000) if x.startswith(':Tags:')]
            if data:
                keywords = [x.strip() for x in data[0][6:].split(' ')]
                for x in keywords:
                    if x:
                        map_keyword2script[x].append(s)

    return map_keyword2script


def printListInColumns(l, ncolumns):
    '''Output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % ncolumns != 0:
        n += 1

    # Build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # Add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # Convert to rows
    rows = list(zip(*columns))

    # Build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for _ in range(ncolumns)])

    # Put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):
    argv = sys.argv

    path = os.path.join(os.path.abspath(os.path.dirname(cgat.__file__)),
                        "tools")

    if len(argv) == 1 or argv[1] in ("--help", "-h"):
        print(globals().get("__doc__", "No documentation available."))

        map_keyword2script = mapKeyword2Script(path)

        if len(argv) <= 2:
            print('cgat tools are grouped by keywords. The following keywords')
            print('are defined:\n')
            print("%s\n" % printListInColumns(list(map_keyword2script.keys()), 3))

        if 'all' in argv[2:]:
            print("The list of all available commands is:\n")
            print("%s\n" % printListInColumns(
                sorted([os.path.basename(x)[:-3]
                        for x in glob.glob(os.path.join(path, "*.py"))]),
                3))
        else:
            for arg in argv[2:]:
                if arg in map_keyword2script:
                    print("Tools matching the keyword '%s':\n" % arg)
                    print('%s\n' % printListInColumns(
                        sorted(map_keyword2script[arg]),
                        3))
        return

    if len(argv) < 2:
        print("Error: No command provided.\nUse --help for more information.", file=sys.stderr)
        sys.exit(1)

    command = argv[1]

    # Replace hyphens with underscores to match Python module naming conventions
    command = re.sub("-", "_", command)

    # Dynamically load the specified module using importlib
    module_path = os.path.join(path, f"{command}.py")

    if not os.path.isfile(module_path):
        print(f"Error: Module '{command}' not found in path '{path}'.", file=sys.stderr)
        sys.exit(1)

    spec = importlib.util.spec_from_file_location(command, module_path)
    if spec is None:
        print(f"Error: Cannot create a module spec for '{command}'.", file=sys.stderr)
        sys.exit(1)

    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        print(f"Error: Failed to load module '{command}': {e}", file=sys.stderr)
        sys.exit(1)

    # Remove 'cgat' from sys.argv to pass the remaining arguments to the module's main function
    if len(sys.argv) > 0:
        del sys.argv[0]

    # Ensure the module has a 'main' function
    if not hasattr(module, 'main'):
        print(f"Error: The module '{command}' does not have a 'main' function.", file=sys.stderr)
        sys.exit(1)

    # Call the module's main function with the remaining arguments
    try:
        module.main(sys.argv)
    except Exception as e:
        print(f"Error: An exception occurred while executing the 'main' function of '{command}': {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
