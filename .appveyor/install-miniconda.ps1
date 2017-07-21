# Sample script to install anaconda under windows
# Authors: Stuart Mumford
# Borrowed from: Olivier Grisel and Kyle Kastner
# License: BSD 3 clause

$MINICONDA_URL = "https://repo.continuum.io/miniconda/"

# We pin the version for conda as it's not the most stable package from
# release to release. Add note here if version is pinned due to a bug upstream.
if (! $env:CONDA_VERSION) {
   $env:CONDA_VERSION = "4.3.6"
}

function DownloadMiniconda ($version, $platform_suffix) {
    $webclient = New-Object System.Net.WebClient
    $filename = "Miniconda3-" + $version + "-Windows-" + $platform_suffix + ".exe"

    $url = $MINICONDA_URL + $filename

    $basedir = $pwd.Path + "\"
    $filepath = $basedir + $filename
    if (Test-Path $filename) {
        Write-Host "Reusing" $filepath
        return $filepath
    }

    # Download and retry up to 3 times in case of network transient errors.
    Write-Host "Downloading" $filename "from" $url
    $retry_attempts = 2
    for($i=0; $i -lt $retry_attempts; $i++){
        try {
            $webclient.DownloadFile($url, $filepath)
            break
        }
        Catch [Exception]{
            Start-Sleep 1
        }
   }
   if (Test-Path $filepath) {
       Write-Host "File saved at" $filepath
   } else {
       # Retry once to get the error message if any at the last try
       $webclient.DownloadFile($url, $filepath)
   }
   return $filepath
}

function InstallMiniconda ($miniconda_version, $architecture, $python_home) {
    Write-Host "Installing miniconda" $miniconda_version "for" $architecture "bit architecture to" $python_home
    if (Test-Path $python_home) {
        Write-Host $python_home "already exists, skipping."
        return $false
    }
    if ($architecture -eq "x86") {
        $platform_suffix = "x86"
    } else {
        $platform_suffix = "x86_64"
    }
    $filepath = DownloadMiniconda $miniconda_version $platform_suffix
    Write-Host "Installing" $filepath "to" $python_home
    $args = "/InstallationType=AllUsers /S /AddToPath=1 /RegisterPython=1 /D=" + $python_home
    Write-Host $filepath $args
    Start-Process -FilePath $filepath -ArgumentList $args -Wait -Passthru
    #Start-Sleep -s 15
    if (Test-Path $python_home) {
        Write-Host "Miniconda $miniconda_version ($architecture) installation complete"
    } else {
        Write-Host "Failed to install Python in $python_home"
        Exit 1
    }
}

# Install miniconda, if no version is given use the latest
if (! $env:MINICONDA_VERSION) {
   $env:MINICONDA_VERSION="latest"
}

InstallMiniconda $env:MINICONDA_VERSION $env:PLATFORM $env:PYTHON

# Set environment variables
$env:PATH = "${env:PYTHON};${env:PYTHON}\Scripts;" + $env:PATH

# Conda config
conda config --set always_yes true
conda config --add channels defaults

if ($env:CONDA_CHANNELS) {
   $CONDA_CHANNELS=$env:CONDA_CHANNELS.split(" ")
   foreach ($CONDA_CHANNEL in $CONDA_CHANNELS) {
           conda config --add channels $CONDA_CHANNEL
   }
   Remove-Variable CONDA_CHANNELS
   rm env:CONDA_CHANNELS
}

# Install the build and runtime dependencies of the project.
conda install -q conda=$env:CONDA_VERSION

if (! $env:CONDA_CHANNEL_PRIORITY) {
   $CONDA_CHANNEL_PRIORITY="false"
} else {
   $CONDA_CHANNEL_PRIORITY=$env:CONDA_CHANNEL_PRIORITY.ToLower()
}

# We need to add this after the update, otherwise the ``channel_priority``
# key may not yet exists
conda config  --set channel_priority $CONDA_CHANNEL_PRIORITY

# Create a conda environment 
conda create -q -n test python=$env:PYTHON_VERSION
activate test

# Set environment variables for environment (activate test doesn't seem to do the trick)
$env:PATH = "${env:PYTHON}\envs\test;${env:PYTHON}\envs\test\Scripts;${env:PYTHON}\envs\test\Library\bin;" + $env:PATH

# Check that we have the expected version of Python
python --version

# CORE DEPENDENCIES: PIP
# Could install other core dependencies here, e.g. pytest
#conda install -q -n test pip pytest
conda install -q -n test pip

# Check whether a specific version of Numpy is required
if ($env:NUMPY_VERSION) {
    if($env:NUMPY_VERSION -match "stable") {
        $NUMPY_OPTION = "numpy=" + $env:LATEST_NUMPY_STABLE
    } elseif($env:NUMPY_VERSION -match "dev") {
        $NUMPY_OPTION = "Cython pip".Split(" ")
    } else {
        $NUMPY_OPTION = "numpy=" + $env:NUMPY_VERSION
    }
    conda install -n test -q $NUMPY_OPTION
} else {
    $NUMPY_OPTION = ""
}

# Install the specified versions of numpy and other dependencies
if ($env:CONDA_DEPENDENCIES) {
    $CONDA_DEPENDENCIES = $env:CONDA_DEPENDENCIES.split(" ")
} else {
    $CONDA_DEPENDENCIES = ""
}

# Check whether the installation is successful, if not abort the build
$output = cmd /c conda install -n test -q $NUMPY_OPTION $CONDA_DEPENDENCIES 2>&1

echo $output
if ($output | select-string UnsatisfiableError, PackageNotFoundError) {
   echo "Installing dependencies with conda was unsuccessful, using pip instead"
   $output = cmd /c pip install $CONDA_DEPENDENCIES 2>&1
   echo $output
   if ($output | select-string UnsatisfiableError, PackageNotFoundError) {
      $host.SetShouldExit(1)
   }
}

# Check whether the developer version of Numpy is required and if yes install it
if ($env:NUMPY_VERSION -match "dev") {
   Invoke-Expression "${env:CMD_IN_ENV} pip install git+https://github.com/numpy/numpy.git#egg=numpy --upgrade --no-deps"
}

# We finally install the dependencies listed in PIP_DEPENDENCIES. 

if ($env:PIP_FLAGS) {
    $PIP_FLAGS = $env:PIP_FLAGS.split(" ")
} else {
    $PIP_FLAGS = ""
}

if ($env:PIP_DEPENDENCIES) {
    $PIP_DEPENDENCIES = $env:PIP_DEPENDENCIES.split(" ")
} else {
    $PIP_DEPENDENCIES = ""
}

if ($env:PIP_DEPENDENCIES) {
    pip install $PIP_DEPENDENCIES $PIP_FLAGS
}