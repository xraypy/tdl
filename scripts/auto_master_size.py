import errno
import os
import subprocess
import sys
import time

def auto_master(spec_dir):
    spec_sizes = {}
    while True: #spec_sizes == {}:
        print 'Updating...'
        spec_files = os.listdir(spec_dir)
        spec_files = [item for item in spec_files if item.endswith('.spc')]
        for item in spec_files:
            full_item = os.path.join(spec_dir, item)
            edit_size = os.path.getsize(full_item)
            if item not in spec_sizes:
                spec_sizes[item] = edit_size
                subprocess.call([sys.executable,
                                 'C:\\apps\\tdl\\scripts\\spectohdf.py',
                                 full_item, '', '-a'])
            else:
                if spec_sizes[item] < edit_size:
                    spec_sizes[item] = edit_size
                    subprocess.call([sys.executable,
                                    'C:\\apps\\tdl\\scripts\\spectohdf.py',
                                     full_item, '', '-a'])
        time.sleep(20)
    
"""
File Locker
Author: Evan Fosmark
http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python
Last modified: 7.16.2012 by Craig Biwer (cbiwer@uchicago.edu)
"""
 
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
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                os.write(self.fd, os.environ['USERNAME'] + '\n')
                os.write(self.fd, os.environ['COMPUTERNAME'] + '\n')
                os.write(self.fd, time.ctime(time.time()))
                break;
            except OSError as e:
                if e.errno != errno.EEXIST and e.errno != errno.EACCES:
                    raise 
                if (time.time() - start_time) >= self.timeout:
                    if e.errno == errno.EEXIST:
                        raise FileLockException("Timeout occured.")
                    else:
                        raise FileLockException("Access denied.")
                time.sleep(self.delay)
        self.is_locked = True
 
 
    def release(self):
        """ Get rid of the lock by deleting the lockfile. 
            When working in a `with` statement, this gets automatically 
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
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


if __name__ == '__main__':
    args = sys.argv[1:]
    auto_master(*args)
