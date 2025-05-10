[![Build package](https://github.com/r-universe/krystianll/actions/workflows/build.yml/badge.svg?branch=master)](https://github.com/r-universe/krystianll/actions/workflows/build.yml)
# mlemur
This project was supported by National Science Centre, Poland (PRELUDIUM Grant 2019/35/N/NZ1/03402).

## Simplest way to install if you have an Intel machine (Windows, Intel-based older Macs)
1. Install newest version of R (currently 4.3.0) in the default location:
- For Intel Macs use the following link: <https://cran.r-project.org/bin/macosx/base/>
- For Windows machines use the following link: <https://cran.r-project.org/bin/windows/base/>

2. Windows only: Install the corresponding version of RTools (e.g., if R version is 4.3.x, then RTools must be 4.3) in the default location from <https://cran.r-project.org/bin/windows/Rtools/>

3. Locate and run R.exe or R.app

4. In the R main window, copy the following command to R console and press Enter:
```{r eval=FALSE, warning=TRUE}
install.packages('mlemur', repos = c('https://krystianll.r-universe.dev', 'https://cloud.r-project.org'))
```

## Instructions (Windows)
1. Install newest version of R (currently 4.3.0) in the default location from <https://cran.r-project.org/bin/windows/base/>

2. Install the corresponding version of RTools (e.g., if R version is 4.3.x, then RTools must be 4.3) in the default location from <https://cran.r-project.org/bin/windows/Rtools/>

3. Proceed to the next part. You can either follow instruction for installing the pre-compiled binary (recommended) or compiling from the source file.

### Installation using a pre-compiled binary package (recommended)
1. Download the latest release of mlemur (zip format) from <https://github.com/krystianll/mlemur/releases/>

2. Locate and run R.exe

3. In the R main window, install required packages: copy the following command to R console and press Enter:
```{r eval=FALSE, warning=TRUE}
install.packages(c("devtools", "Rcpp", "reactable", "readxl", "writexl", "shiny",
                   "shinyFeedback", "shinyWidgets", "shinyjs", "boot", "BH",
                   "hypergeo", "Ryacas"), dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

4. In the R main window, select "Packages" in the menu bar at the top of the window and then "Install package(s) from local zip files"

5. In the new window, locate the downloaded .zip file and click Open

### Installation using source files (recommended; requires compilation)
1. Locate and run R.exe

2. Install `devtools` package: copy the following command to R console and press Enter:
```{r eval=FALSE}
install.packages("devtools", dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

3. Install mlemur with all dependencies using the following command:
```{r eval=FALSE}
devtools::install_github("krystianll/mlemur", ref = "master",
                         dependencies = TRUE)
```

### Running mlemur
Mlemur can be run in both text and graphical mode. The graphical mode is initiated in a browser window.

To initiate the graphical mode, use the following command in R:
```{r eval=FALSE}
mlemur::mlemur()
```
Alternatively, you can use the following clickable script: <https://github.com/krystianll/mlemur/blob/master/run_mlemur.bat>
1. Copy the text within the file (everything between and including lines 1 and 35).

2. Paste the text into a new empty Notepad document.

3. Save the file as run_mlemur.bat. If asked to confirm the the change in file extension, press OK. (Make sure the filename is not run_mlemur.bat.txt!)

## Instructions (macOS)
1. Install newest version of R (currently 4.3.0) in the default location.
- For Intel Macs use the following link: <https://cran.r-project.org/bin/macosx/base/>
- For Apple Silicon (M1, M2) Macs use the following link: <https://cran.r-project.org/bin/macosx/big-sur-arm64/base/>
- If you're not sure which version to choose, click Apple logo in the top-left corner of the screen, choose About This Mac and inspect the Chip model

2. Don't open R yet. Open Terminal.app. To do this, press and hold Command ⌘ key and then press Space key. Type Terminal and press Enter.

3. In the Terminal window, copy the following command and press Enter:
```{r eval=FALSE}
xcode-select --install
```
Then follow the instructions in the new window.

4. Proceed to the next part. You can either follow instruction for installing the pre-compiled binary or compiling from the source file.

### Installation using a pre-compiled binary package (recommended)
1. Locate and run R.app

2. Download the latest release of mlemur (.tgz format)corresponding to the verion of R you have installed from the following link: <https://github.com/krystianll/mlemur/releases/>.
- For Intel Macs use the file with _x64 in the filename
- For Apple Silicon (M1, M2) Macs use the file with _arm in the filename
- If you're not sure which version to choose, inspect the text in R main window:
- Platform: aarch64-apple-darwinXX (64-bit) (with XX some number) means you should download the arm version
- Platform: x86_64-apple-darwinXX (64-bit) (with XX some number) means you should download the intel version

**Make sure that the file has .tgz extension! If it is saved as .tar, please change the extension by hand.**

3. In the R main window, install required packages: copy the following command to R console and press Enter:
```{r eval=FALSE}
install.packages(c("devtools", "Rcpp", "reactable", "readxl", "writexl", "shiny",
                   "shinyFeedback", "shinyWidgets", "shinyjs", "boot", "BH",
                   "hypergeo", "Ryacas"), dependencies = TRUE)
```
You might be asked to confirm by clicking "Yes" or to select server from which the packages will be downloaded - select the first one

4. In the R main window, select "Packages & Data" in the menu bar at the top of the window and then "Package Installer"

5. In the top-left menu, in the "Packages repository" section change CRAN (binaries) to Local Binary Package and then click Install in the bottom-right part of the window

6. Locate the mlemur .tgz file and press Open

Alternative 4-6. If the above doesn't work, in the R main window, select "Misc" in the menu bar at the top of the window and then "Change Working Directory…". Locate the folder containing the .tgz file and press Open. Then in the R console use the following command (replace the name_of_the_package with a valid filename) and press Enter:
```{r eval=FALSE}
install.packages("name_of_the_package.tgz")
```

### Installation using source files (requires compilation)
1. Open Terminal.app. To do this, press and hold Command ⌘ key and then press Space key. Type Terminal and press Enter.

2. In the Terminal window, copy the following command and press Enter to install homebrew:
```{r eval=FALSE}
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
You might be prompted to type your password and press Enter. Press Enter again when asked to accept the creation of new folders. Installation might take several minutes.

3. In the Terminal window, copy the following command and press Enter to install boost libraries:
```{r eval=FALSE}
brew install boost
```
4. In the Terminal window, copy the following command and press Enter to install mpfr:
```{r eval=FALSE}
brew install mpfr
```
5. Locate and run R.app
6. Install `devtools` package: copy the following command to R console and press Enter:
```{r eval=FALSE}
install.packages("devtools", dependencies = TRUE)
```
You might be asked to confirm by clicking 'Yes' or to select server from which the packages will be downloaded - select the first one

7. Install mlemur with all dependencies using the following command:
```{r eval=FALSE}
devtools::install_github("krystianll/mlemur", ref = "master",
                         dependencies = TRUE)
```

### Running mlemur
Mlemur can be run in both text and graphical mode. The graphical mode is initiated in a browser window.

To initiate the graphical mode, use the following command in R:
```{r eval=FALSE}
mlemur::mlemur()
```
Alternatively, you can use the following clickable script: <https://github.com/krystianll/mlemur/blob/master/run_mlemur.command>
1. Copy the text within the file (lines 1 and 2).

2. Paste the text into a new empty TextEdit document.

3. In the Menu bar, click Format and choose "Make Plain Text".

4. Save the file as run_mlemur.command. If asked to confirm the the change in file extension, press OK.

5. Open the Terminal window and copy the following command:
```{r eval=FALSE}
chmod u+x 
```
Make sure there is a space after "x". Drag and drop the run_mlemur.command file to Terminal window and then press Enter.

6. Open the run_mlemur.command file. macOS might still try to block the execution of the file. If that is the case, click on the Apple logo in the top-left corner of the screen, select System Preferences, then Security & Privacy. Click "Open Anyway" and then confirm by clicking "Open". The file will work from now on.
