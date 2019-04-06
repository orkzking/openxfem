# openxfem
Automatically exported from [code.google.com/p/openxfem](https://code.google.com/archive/p/openxfem)

## What is this Project for?
This is an opensource code I found on google code several years ago.
As the name let one guess it is a eXtended Finite Element Method (**XFEM**) code.
As far as I looked through the source code, it doesn't run in parallel and has got some runtime issues (only compile in debug mode).
At the moment I use this code to test the functionallities of Git & GitHub.
E.g. working with a gh-pages branch to easily provide the documentation of this code online [here](https://orkzking/github.io/openxfem).
Also a Gitpitch presentation is in progress [here](https://gitpitch.com/orkzking/openxfem/master).
This is planed to be a step-by-step quickstart guide.
## How to configure and install?
To get this code runnig I created branches for different OS and/or IDEs.
The nearly original code in the master branch seems to compile under Windows with MinGW (I never tested this).
The following branches are tested to compile with the namegiving IDE on a certain OS and containing a project of this IDE.
 - **visualstudio2010**: Is tested for Win 8.1 and visual studio 2010 ultimate with the git plugin [Git Source Control Provider](https://marketplace.visualstudio.com/items?itemName=yysun.GitSourceControlProvider)
 - **eclipse**: Is tested for Debian and Eclipse Oxygen with the CDT-plugin and gcc toolchain; For now the changes are quick & dirty to get it run

 The idea is to "hopefully" just checkout the branch of choice, open the IDE and press the build button. :joy:
 More branches are welcome e.g. one is planed for Termux on Android.

## The License of the Project
[GNU GPL v3](https://gnu.org/licenses/gpl.html)
