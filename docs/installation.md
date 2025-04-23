# Installation Guide

This guide explains how to install the VCFX toolkit on your system.

## Prerequisites

Before installing VCFX, ensure you have the following prerequisites:

- A Unix-like operating system (Linux, macOS)
- C++ compiler (GCC 5.0+ or Clang 3.8+)
- CMake (version 3.10 or higher)
- Git (for cloning the repository)

## Installation Methods

### Method 1: Using Bioconda (Recommended)

VCFX is available in the Bioconda channel, providing an easy installation method:

```bash
# Set up Bioconda channels if you haven't already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install VCFX
conda install vcfx
```

This method automatically handles all dependencies and provides a ready-to-use installation of VCFX.

### Method 2: Building from Source

This method ensures you have the latest version of VCFX.

1. Clone the VCFX repository:

   ```bash
   git clone https://github.com/ieeta-pt/VCFX.git
   cd VCFX
   ```

2. Create a build directory and build the tools:

   ```bash
   mkdir -p build
   cd build
   cmake ..
   make
   ```

3. (Optional) Install the tools to your system:
   ```bash
   sudo make install
   ```

After installation, you should be able to run VCFX tools from your terminal.

### Method 3: Building Individual Tools

If you only need specific tools, you can build them individually:

```bash
cd VCFX/build
make VCFX_tool_name
```

Replace `VCFX_tool_name` with the actual tool name (e.g., `VCFX_allele_freq_calc`).

### Method 4: Using Docker (Coming Soon)

We are working on providing Docker images for VCFX to simplify the installation process.

## Verifying Installation

To verify that VCFX tools are installed correctly, you can run:

```bash
VCFX_validator --help
```

You should see the help message for the VCFX_validator tool.

## Adding VCFX to Your PATH

If you didn't use `make install` or want to run tools directly from the build directory, you can add them to your PATH:

```bash
# Add to ~/.bashrc or ~/.zshrc
export PATH="$PATH:/path/to/VCFX/build/src"
```

Replace `/path/to/VCFX` with the absolute path to your VCFX repository.

## Troubleshooting

### Missing Dependencies

If CMake complains about missing dependencies, you may need to install additional libraries:

**Ubuntu/Debian:**
```bash
sudo apt-get install build-essential cmake git
```

**macOS (with Homebrew):**
```bash
brew install cmake
```

### Compilation Errors

If you encounter compilation errors, ensure you're using a compatible compiler version and have all necessary dependencies installed.

### Path Issues

If you get "command not found" errors when trying to run VCFX tools, check that:

1. The tools were built successfully
2. The tools are in your PATH
3. The executable files have execute permissions

## Getting Help

If you encounter any issues during installation, please:

1. Check the [GitHub Issues](https://github.com/ieeta-pt/VCFX/issues) to see if your problem has been reported
2. Open a new issue if necessary, including details about your system and the exact error messages

## Next Steps

Now that you have VCFX installed, you can:

- Follow the [Quick Start Guide](quickstart.md) to learn how to use VCFX
- Explore the [Tool Documentation](tools_overview.md) to learn about specific tools 