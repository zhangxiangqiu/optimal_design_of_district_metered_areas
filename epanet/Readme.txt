Contents of TOOLKIT.ZIP
=======================

EPANET is a program that analyzes the hydraulic and water quality behavior of water
distribution systems. The EPANET Programmer's Toolkit is a dynamic link library (DLL) of
functions that allows developers to customize EPANET's computational engine for their 
own specific needs. The functions can be incorporated into 32-bit Windows applications 
written in C/C++, Delphi Pascal, Visual Basic, or any other language that can call
functions within a Windows DLL. The Toolkit DLL file is named EPANET2.DLL and is
distributed with EPANET. The Toolkit comes with several different header files, function 
definition files, and .lib files that simplify the task of interfacing it with C/C++, 
Delphi, and Visual Basic code.

This archive contains the following library files:
  EPANET2.DLL   -- the Toolkit function library
  EPANET2.LIB   -- LIB file for using the Toolkit with Microsoft Visual C++

epanet2.lib must be linked in to any Toolkit application compiled for Windows using MS Visual C++.

The Toolkit provides several header files that are needed to develop C/C++ applications:

epanet2.h contains declarations of the single-threaded version of the Toolkit (the ENxxx named functions).
epanet2_2.h contains declarations of the multi-threaded version of the Toolkit (the EN_xxx named functions).
epanet2_enums.h contains definitions of the symbolic constants used by the Toolkit.

Developers need to issue an include directive for either epanet2.h or epanet2_2.h in their C/C++ code depending on whether they are building 
a single-threaded or multi-threaded application. There is no need to explicitly include epanet2_enums.h as it is automatically included by 
both of the other header files.

Several additional function declaration files that provide bindings for other programming languages are included with the Toolkit package:

epanet2.bas for Visual Basic for Applications and Visual Basic 6
epanet2.vb for Visual Basic .NET
epanet2.pas for Delphi Pascal, Free Pascal or Lazarus.
These bindings only support the single-threaded version of the Toolkit.

More infomation about the Toolkit and the API reference can be accessed at: http://wateranalytics.org/EPANET/


