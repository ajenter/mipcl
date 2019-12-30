
Welcome to MIPCL
================

The MIPCL installation bundle includes the MIPCL itself (the library, headers, documentation, and examples),
a modeling tool called MIPshell,
and MIPCL-PY, which is a collection of modules that allow us using MIPCL in Python programs.
All information about MIPCL, MIPshell, and MIPCL-PY is available on www.mipcl-cpp.appspot.com

Installing MIPCL
----------------

1. Download a MIPCL archive from www.mipcl-cpp.appspot.com.
   The site provides download links for all supported development platforms.

2. Unpack the MIPCL archive (tar.gz package) into your installation dirertory
  (for example, /usr/local, /opt, or any subdirectory of your home directory).

3. Enter the top-level MIPCL-directory "mipcl-<version>" (for example, <version>=1.4.0):
   cd mipcl-<version>

4. Run the following script with superuser privileges:
   ./install_libs
    Alternatively, you can
       a) create symbolic links:
        cd lib
        ln -s libmipcl.so.<version> libmipcl.so
        ln -s libmipshell.so.<version> libmipshell.so
        (Optionally, if you will use JVmipcl)
        ln -s libJVmipcl.so.<version> libJVmipcl.so
        (Optionally, if you will use MIPCL-PY)
        cd ../mipcl_py/mipshell
        If you use Python 2,
        ln -s mipcl-py2.so.<version> mipcl.so
        or if you use Python 3,
        ln -s mipcl-py3.so.<version> mipcl.so   
     and b) modify the LD_LIBRARY_PATH environment variable
        to include the path to the MIPCL libraries.
        If you use bash, add the following line to the end of .bash_profile:
        LD_LIBRARY_PATH=$MIPCLDIR/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
 
5. Define a new environment variable MIPCLDIR that points to <install_dir>/mipcl-<version>,
   where <install_dir> refers to the directory that contains the "mipcl-<version>" directory,
   and (optionally) alter the $PATH variable to include the path to the MIPCL binaries.
   For example, if you use the bash shell, add the following to the end of your .bash_profile:
   MIPCLDIR="<install_dir>/mipcl-<version>"; export MIPCLDIR
   PATH=$MIPCLDIR/bin:$PATH; export PATH

6. (Optional) If you have installed MIPCL into a directory not avalable for writing,
    copy (or, if you are the only developer using MIPCL, move)
    the `examples` directory to any place in your home directory.

REMARK:
   ~/.bash_profile or ~/.profile files are not executed
   when logging in into a Wayland session.
   It seams that there is no satisfactory way of setting environment variables
   when logging in. A few remedies are still possible. For example,
   you can move a file with the following content

MIPCLDIR=<install_dir>
PYTHONPATH=$MIPCLDIR
PATH=$PATH:$MIPCLDIR/bin
LD_LIBRARY_PATH=$MIPCLDIR/lib

   into /usr/share/gdm/env.d directory.
   To do this, you need administrator privileges.

   
   When you run your programs from the console,
   you can just move the content of ~/.bash_profile into ~/.bashrc,
   and then, anytime a new shell session starts, ~/.bashrc is executed
   and all you environment  variables are set.
