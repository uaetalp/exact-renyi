"""Portable resource checks for standalone, single-process experiment workers."""
import os
import resource
import sys
import time

_deadline = None


def configure(seconds=180):
    global _deadline
    _deadline = time.monotonic() + seconds
    cpu_limit = max(1, int(seconds))
    resource.setrlimit(resource.RLIMIT_CPU, (cpu_limit, cpu_limit))
    if sys.platform.startswith('linux'):
        resource.setrlimit(resource.RLIMIT_AS, (2**31, 2**31))
        if hasattr(os, 'sched_getaffinity'):
            os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})


def check_resources():
    if _deadline is not None and time.monotonic() > _deadline:
        raise TimeoutError('worker time limit exceeded')
    if sys.platform.startswith('linux'):
        if len(os.listdir('/proc/self/task')) != 1:
            raise RuntimeError('experiments require one worker thread')


def peak_rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == 'darwin' else value * 1024
