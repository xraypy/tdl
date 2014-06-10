"""
File Locker
Author: Evan Fosmark
http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python
Last modified: 7.16.2012 by Craig Biwer (cbiwer@uchicago.edu)
"""


import os
import time
import errno

class FileLockException(Exception):
    pass

class FileLock(object):
    """ A file locking mechanism that has context-manager support so
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
    """

    def __init__(self, file_name, timeout=10, delay=.05):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        #self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.lockfile = file_name + '.lock'
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay


    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws
            an exception.
        """
        start_time = time.time()
        lockfile = os.path.abspath(self.lockfile)
        lockdir  = os.path.dirname(lockfile)
        if not os.access(lockdir, os.W_OK):
            raise FileLockException("Cannot write to %s." % lockdir)

        dat = [os.environ.get('USERNAME', ''),
               os.environ.get('COMPUTERNAME', '')]
        while (time.time() - start_time) < self.timeout:
            if os.path.exists(lockfile):
                time.sleep(self.delay)
                continue
            dat.extend([time.ctime(), ''])
            with open(self.lockfile, 'w') as fh:
                fh.write('\n'.join(dat))
            self.is_locked = True
            return
        raise FileLockException("Timeout occured waiting for lockfile.")

    def release(self):
        """ Get rid of the lock by deleting the lockfile.
            When working in a `with` statement, this gets automatically
            called at the end.
        """
        if self.is_locked:
            os.unlink(self.lockfile)
            self.is_locked = False


    def __enter__(self):
        """ Activated when used in the with statement.
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self


    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()


    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()
