'''
bed.plot.py - create genomic snapshots using the IGV Viewer
===========================================================

:Tags: Python

Purpose
-------

Create genomic plots in a set of intervals using
the IGV snapshot mechanism.

The script can use a running instance of IGV identified
by host and port. Alternatively, it can start IGV and load
a pre-built session.

Usage
-----

Example::

   python bed2plot.py < in.bed

Type::

   python script_template.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import socket
import pysam


import cgatcore.experiment as E


class IGV(object):
    """based on IGV.py by Brent Petersen, see here:
    https://github.com/brentp/bio-playground/blob/master/igv/igv.py

    (MIT licenced)
    """

    _socket = None
    _path = None

    def __init__(self, host='127.0.0.1', port=60151, snapshot_dir='/tmp/igv'):
        self.host = host
        self.port = port
        self.commands = []
        self.connect()
        self.set_path(snapshot_dir)

    def connect(self):
        if self._socket:
            self._socket.close()
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect((self.host, self.port))

    def go(self, position):
        return self.send('goto ' + position)
    goto = go

    def genome(self, name):
        return self.send('genome ' + name)

    def load(self, url):
        return self.send('load ' + url)

    def region(self, contig, start, end):
        return self.send(' '.join(map(str, ['region', contig, start, end])))

    def sort(self, option='base'):
        """
        options is one of: base, position, strand, quality, sample, and
        readGroup.
        """
        assert option in ("base", "position", "strand", "quality", "sample",
                          "readGroup")
        return self.send('sort ' + option)

    def set_path(self, snapshot_dir):
        if snapshot_dir == self._path:
            return
        if not os.path.exists(snapshot_dir):
            os.makedirs(snapshot_dir)

        self.send('snapshotDirectory %s' % snapshot_dir)
        self._path = snapshot_dir

    def expand(self, track=''):
        self.send('expand %s' % track)

    def collapse(self, track=''):
        self.send('collapse %s' % track)

    def clear(self):
        self.send('clear')

    def send(self, cmd):
        # socket in Python2 oprates with strings
        if sys.version_info.major == 2:
            self._socket.send(cmd + '\n')
            return self._socket.recv(4096).rstrip('\n')
        # while socket in Python3 requires bytes
        else:
            self.commands.append(cmd)
            cmd = cmd + '\n'
            self._socket.send(cmd.encode('utf-8'))
            return self._socket.recv(4096).decode('utf-8').rstrip('\n')

    def save(self, path=None):
        if path is not None:
            # igv assumes the path is just a single filename, but
            # we can set the snapshot dir. then just use the filename.
            dirname = os.path.dirname(path)
            if dirname:
                self.set_path(dirname)
            return self.send('snapshot ' + os.path.basename(path))
        else:
            return self.send('snapshot')
    snapshot = save


def main(argv=sys.argv):

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    parser.add_option("-s", "--session", dest="session",
                      type="string",
                      help="load session before creating plots "
                      "[%default]")

    parser.add_option("-d", "--snapshot-dir", dest="snapshotdir",
                      type="string",
                      help="directory to save snapshots in [%default]")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("png", "eps", "svg"),
                      help="output file format [%default]")

    parser.add_option("-o", "--host", dest="host", type="string",
                      help="host that IGV is running on [%default]")

    parser.add_option("-p", "--port", dest="port", type="int",
                      help="port that IGV listens at [%default]")

    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend each interval by a number of bases "
                      "[%default]")

    parser.add_option("-x", "--expand", dest="expand", type="float",
                      help="expand each region by a certain factor "
                      "[%default]")

    parser.add_option("--session-only", dest="session_only",
                      action="store_true",
                      help="plot session after opening, "
                      "ignore intervals "
                      "[%default]")

    parser.add_option("-n", "--name", dest="name", type="choice",
                      choices=("bed-name", "increment"),
                      help="name to use for snapshot "
                      "[%default]")

    parser.set_defaults(
        command="igv.sh",
        host='127.0.0.1',
        port=61111,
        snapshotdir=os.getcwd(),
        extend=0,
        format="png",
        expand=1.0,
        session=None,
        session_only=False,
        keep_open=False,
        name="bed-name",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv, add_output_options=True)

    igv_process = None
    if options.new_instance:
        E.info("starting new IGV process")
        igv_process = IGV.startIGV(command=options.command,
                                   port=options.port)
        E.info("new IGV process started")

    E.info("connection to process on %s:%s" % (options.host, options.port))
    E.info("saving images in %s" % options.snapshotdir)
    igv = IGV(host=options.host,
              port=options.port,
              snapshot_dir=os.path.abspath(options.snapshotdir))

    if options.session:
        E.info('loading session from %s' % options.session)
        igv.load(options.session)
        E.info('loaded session')

    if options.session_only:
        E.info('plotting session only ignoring any intervals')
        fn = "%s.%s" % (os.path.basename(options.session), options.format)
        E.info("writing snapshot to '%s'" %
               os.path.join(options.snapshotdir, fn))
        igv.save(fn)

    else:
        c = E.Counter()
        for bed in pysam.tabix_iterator(options.stdin,
                                        parser=pysam.asBed()):

            c.input += 1

            # IGV can not deal with white-space in filenames
            if options.name == "bed-name":
                name = re.sub("\s", "_", bed.name)
            elif options.name == "increment":
                name = str(c.input)

            E.info("going to %s:%i-%i for %s" %
                   (bed.contig, bed.start, bed.end, name))

            start, end = bed.start, bed.end
            extend = options.extend
            if options.expand:
                d = end - start
                extend = max(extend, (options.expand * d - d) // 2)

            start -= extend
            end += extend

            igv.go("%s:%i-%i" % (bed.contig, start, end))

            fn = E.get_output_file("%s.%s" % (name, options.format))
            E.info("writing snapshot to '%s'" % fn)
            igv.save(fn)

            c.snapshots += 1

        E.info(c)

    if igv_process is not None and not options.keep_open:
        E.info('shutting down IGV')
        igv_process.send_signal(signal.SIGKILL)
        
    E.stop()

if __name__ == "__main__":
    sys.exit(main())
