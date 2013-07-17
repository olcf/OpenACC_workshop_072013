A simple script that gathers information about the attached GPU.

See the header at the top of the file practical0.wlm for details.

JOB SUBMISSION:
===============

bash submit.bash

This will:

* load the correct modules using script ../Tools/XK_setup.bash
* create a directory on the lustre filesystem (and report its name)
* write and submit a WLM jobscript in this directory.
  -  Type "qstat -u $USER" or "squeue -u $USER" to check progress and completion.

The directory name is reported. You can then cd to this directory and look at
the output. The output will be in a file with filename ending ".log".

EOF
