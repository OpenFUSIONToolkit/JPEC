# JPEC

This is a work in progress juliafication of the Generalized Perturbed Equilibrium Code suite https://github.com/PrincetonUniversity/GPEC.

## Documentation

The WIP documentation can be found [here](https://openfusiontoolkit.github.io/JPEC/dev/).
## Developer Notes

### Commit messages

To assist with the process of compiling release notes, please confirm to the commit message format:
```
CODE - TAG - Detailed message
```
where CODE is EQUIL, DCON, VAC, etc. and TAGs are short descriptors of the type of commit. Example tags might be WIP (work in progress), MINOR, IMPROVEMENT, BUG FIX, NEW FEATURE, etc. Look through old commits for examples of what tags are commonly used. 

Again, this is currently used for the by-hand compilation of release notes. The tags thus need to be human readable but are not strictly inforced to be within some limited set. The objective is to allow a lead developer to skim through commits and pick out only the key new features / bug fixes / improvements to note in the release while not having to read all the work-in-progress or minor changes.

### Code Formatting

To retain consistency in code formatting, we request that all PRs include code formatting cleanup before merging. We will make use of [JuliaFormatter.jl] (https://domluna.github.io/JuliaFormatter.jl/stable/), which can be used in both Vim and VSCode. For Vim installation, follow the instructions [here] (https://github.com/kdheepak/JuliaFormatter.vim). In VSCode, the [Julia] (https://www.julia-vscode.org/docs/dev/userguide/formatter/) extension has built-in capabiltiies for JuliaFormatter.jl. The [folder formatter] (https://marketplace.visualstudio.com/items?itemName=MichaKaleta.folderformatter&ssr=false#overview) extension allows formatting of all files within a given folder, which can be useful for larger commits. 

The JuliaFormatter.toml file contains the default settings we will use; this should not be changed unless reformatting the entire repository. 

### Using github

Please see [this link](https://docs.google.com/document/d/1XAOTz1IV8ErZAAk-iSuEuddNOLB5XcoVZsAbPKRUUuA/edit?usp=sharing) for a guide on how to discuss your code development using github. 
