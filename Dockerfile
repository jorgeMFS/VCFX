FROM ubuntu:22.04 AS builder

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a working directory
WORKDIR /app

# Copy the source code
COPY . .

# Build VCFX
RUN mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j$(nproc)

# Create a smaller runtime image
FROM ubuntu:22.04

# Add metadata labels
LABEL org.opencontainers.image.title="VCFX"
LABEL org.opencontainers.image.description="A comprehensive VCF manipulation toolkit"
LABEL org.opencontainers.image.url="https://github.com/ieeta-pt/VCFX"
LABEL org.opencontainers.image.documentation="https://ieeta-pt.github.io/VCFX/"
LABEL org.opencontainers.image.source="https://github.com/ieeta-pt/VCFX"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.citation="Silva, J.M., Oliveira, J.L. (2025). VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis. IWBBIO 2025."

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libz1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the built executables from builder stage
COPY --from=builder /app/build/src /usr/local/bin/

# Create a directory for data
WORKDIR /data

# Add the helper scripts
COPY add_vcfx_tools_to_path.sh /usr/local/bin/
COPY docker_entrypoint.sh /usr/local/bin/

# Make them executable
RUN chmod +x /usr/local/bin/add_vcfx_tools_to_path.sh /usr/local/bin/docker_entrypoint.sh

# Use a custom entrypoint that sets up PATH for the tools
ENTRYPOINT ["/usr/local/bin/docker_entrypoint.sh"]

# Default command shows available tools
CMD ["bash", "-c", "echo 'VCFX Toolkit is ready. Run any VCFX tool by name, for example:' && ls -1 /usr/local/bin/VCFX_* | xargs -n1 basename"]
