# Setting up JPEC 

- [Setting up JPEC](#setting-up-jpec)
  - [On Windows via WSL (Ubuntu)](#on-windows-via-wsl-ubuntu)
  - [On macOS](#on-macos)

## On Windows via WSL (Ubuntu)
1. Install WSL and Ubuntu
   If you don't already have WSL and Ubuntu installed, set this up. [This page](https://learn.microsoft.com/en-us/windows/wsl/install) gives detailed instructions on how to complete the installation. In the Windows Powershell,
   
   1. Make sure WSL is installed:
        ```PowerShell
        wsl --install
        ```
    2. Set Ubuntu as your default WSL distro:
        ```PowerShell
        wsl --set-default Ubuntu
        ```
    3. Launch Ubuntu and update:
        ```PowerShell
        sudo apt update && sudo apt upgrade -y
        ```
2. Install build tools in WSL
    ```Powershell
    sudo apt install build-essential gfortran cmake -y
    ```

    `build-essential` → GCC, make

    `gfortran` → Fortran compiler

    `cmake` → sometimes needed by dependencies

3. Install Julia in WSL
   
    1. Download the latest Linux tarball from the official site [Julia downloads](https://julialang.org/downloads/). It will look like
        ``` shell 
        wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-1.11.3-linux-x86_64.tar.gz
        ```

        ☆ Replace the URL with the latest stable version.

    2. Extract and move it to /opt (or any path):

        ``` shell
        tar -xvzf julia-1.11.3-linux-x86_64.tar.gz
        sudo mv julia-1.11.3 /opt/
        ```
        
        ☆ Ensure these commands match the tarball you installed. These commands match the above tarball and might need to be modified for you installation.
    
    3. Add Julia to PATH:

        ``` shell 
        echo 'export PATH=/opt/julia-1.11.3/bin:$PATH' >> ~/.bashrc
        source ~/.bashrc
        ```

    4. Test it is properly installed

        ``` shell
        julia --version 
        ```

4. Install Python/Jupyter in WSL 
   
   This step is only really required if you want to run the `.ipynb` test notebooks.
   1. To install Python3 and Jupyter notebooks, use these commands
        ``` shell
        sudo apt install python3-pip python3-venv -y
        python3 -m pip install --user jupyter jupyterlab notebook ipykernel
        ```
    ⚠ Important: Add local Python scripts to PATH:

        ``` shell
        echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
        source ~/.bashrc
        ```
    2. Verify it is properly installed

        ``` shell
        jupyter --version
        ```

5. Clone or move JPEC into WSL home folder.
It’s recommended to avoid OneDrive for builds (slower, permission issues):
    ``` shell
    cp -r /mnt/c/Users/<Username>/<Path>/JPEC ~/JPEC
    cd ~/JPEC
    ```

    Where `<Username>/<Path>` is replaced with the path to the JPEC file you cloned to your local machine.

    **OR**

    Alternatively you can clone it from GitHub directly to your virtual machine.

    ``` shell
    git clone https://github.com/OpenFUSIONToolkit/JPEC.git
    cd JPEC
    ```

    Technically, you can use the don't have to copy it from your local installation and you can just `cd mnt/c/Users/<Username>/<Path>/JPEC`, but this in not recommended.

6. Build Fortran dependencies (libspline.so)
    1. Go to the spline source folder:
        ``` shell
        cd ~/JPEC/src/Splines/fortran
        ```

    2. Edit Makefile for Linux/WSL compatibility. SInce many of the current developers use macOS systems, the Makefile is set up for macOS. In the makefile (`JPEC/src/Splines/fortran/makefile`), set
        ``` Makefile
        LIBS ?=
        ```

        instead of 
        ``` Makefile
        LIBS ?=
        ```

        which only works for macOS. Additionally `.dylib` is a file type that is not recognized on Linux, instead Linux uses `.so` files. So, add 
        ``` Makefile
        LIBSUFFIX ?= .so
        ```

        to the makefile in place of 
        ``` Makefile
        LIBSUFFIX ?= .dylib
        ```
        the Mac version.

        ⚠ Important: Be careful **not to commit** changes to the makefile as most developers use macOS.


    3. Back in the shell, clean previous builds using
        ``` shell
        make clean
        ```

    4. Build
        ``` shell
        make
        ```

    5. Verify the library exists using
        ``` shell
        ls ../../../deps/libspline.so
        ```

   6. Export library path so Julia can find it
        ``` shell
        export LD_LIBRARY_PATH=~/JPEC/deps:$LD_LIBRARY_PATH
        ```

        Optional: add to `~/.bashrc` for persistence.

7. Install the Julia packages for JPEC
    1. Launch Julia:
        ``` shell
        julia
        ```
    2. In Julia REPL:
        ``` julia
        using Pkg
        Pkg.instantiate()       # install recorded dependencies
        Pkg.add("Preferences")  # install missing dependency if needed
        Pkg.build("IJulia")     # rebuild kernel
        Pkg.precompile()        # precompile all packages - probably unnecessary 
        ```

8. At this point, you should be able to run the code, open a `.ipynb` notebook, or connect VS Code to your WSL session.
    1. To open a .ipynb notebook
        1. Launch Jupyter from WSL, make sure you have exited Julia using the `exit()` command and then type in the shell
        ```shell
        jupyter notebook --no-browser
        ```
        It will print a URL with a token.

        2. Copy the URL into your Windows browser **OR** open the notebook in VS Code using the **Remote - WSL** extension.

    2. (Optionally) Integrate WSL with VS Code
        1. Install **Remote - WSL** extension in VS Code.
        2. Open VS Code → Connect To → Connect to WSL.
        3. Click Open Folder and then navigate to the JPEC folder on your VM. Open your JPEC folder from WSL: ~/JPEC.
   
        If this is not working, you can launch vscode from the WSL shell you have using the command `code .`
        
        4. Open a terminal inside VS Code — it will automatically use WSL/Ubuntu.
        5. You can now run:
            ``` shell
            make        # rebuild libspline.so if needed
            julia       # run scripts
            jupyter notebook --no-browser
            ```

        6. VS Code also lets you open `.ipynb` notebooks in the WSL environment using the Jupyter extension. Click the "Select Kernel" button in the top right hand of the `.ipynb` file and select the Julia kernel installed in WSL.All dependencies (libspline.so, Julia packages) are accessible.
    3.  Run JPEC
        1. Make sure you are in WSL terminal, with `LD_LIBRARY_PATH` set to include deps.
        2. Launch Julia and run your scripts as usual:
            ``` shell
            include("path/to/jpec_script.jl")
            ```

## On macOS
