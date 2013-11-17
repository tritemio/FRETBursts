To compile the code into C run (requires cython):

    python setup.py build_ext --inplace

Windows 7 x64 requirements
==========================
In windows you need to install a 64-bit C++ compiler. The easiest way is to 
install the MS Windows SDK v7.0 (GRMSDKX_EN_DVD.iso). In principle should not 
be necessary but I installed also Visual Studio 2008 Express SP1
(VS2008ExpressWithSP1ENUX1504728.iso).

After installation and rebooting open the SDK CMD shell (Programs->
Microsoft Windows SDK v7.0->CMD shell) and type:

> set DISTUTILS_USE_SDK=1
> setenv /x64 /release

(now the shell turns green) to configure the  compiler for 64bit output. 

Now, if you have installed a python distribution that set the system-wide %path 
enviroment variable (like Anaconda) then you can compile right away. 

Just cd inside the folder where the pyx file is and run:

    python setup.py build_ext --inplace

Setting the PATH variable
-------------------------

If the %path variable is not configured by your python distribution (like with 
WinPython) you can set the %path manually with "set" command.

For example, with WinPython2.7.3.3, type:

>set path=C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\lib\site-packages\numpy\core;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\Lib\site-packages\PyQt4;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\DLLs;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\Scripts;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\..\tools;C:\Data\Antonio\software\WinPython-64bit-2.7.3.3\python-2.7.3.amd64\..\tools\gnuwin32\bin;C:\Python27\Lib\site-packages\PyQt4;%path%

To obtain this command I just launched "WinPython Command Prompt.exe" typed:

> echo %path%

from the ouput I copied the first part that refers to python.



