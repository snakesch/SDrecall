#! /usr/bin/env python3

import os
import subprocess
import logging
import tempfile



def executeCmd(cmd, logger = logging.getLogger('SDrecall')) -> None:

    proc = subprocess.run(cmd, shell=True, capture_output=True)

    logger.debug(f"Command run:{cmd}")
    logger.debug(f"Return code: {proc.returncode}")

    err_msg = proc.stderr.decode()
    if err_msg != "":
        logger.debug(proc.stderr.decode())

    return proc.stdout.decode()



def prepare_tmp_file(tmp_dir="/paedyl01/disk1/yangyxt/test_tmp", **kwargs):
    try:
        os.mkdir(tmp_dir)
    except FileExistsError:
        pass

    return tempfile.NamedTemporaryFile(dir = tmp_dir, delete = False, **kwargs)