# External tools plus the locked Python dependencies. The repo is installed editable at /app/src and that
# directory is emptied: tasks bind the checkout's src there, so a `git pull` changes the code and only a uv.lock
# change needs a rebuild (and a missing bind fails at import rather than running stale code).
# Ubuntu 24.04 because skiver's release binary is built on ubuntu-latest (glibc 2.39).
FROM ubuntu:24.04
ARG SKIVER=v0.3.2
ARG MINIMAP2=2.31
ARG MINIBWA=v0.7
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl bzip2 git make gcc libc6-dev \
        zlib1g-dev libgomp1 procps \
    && curl -fsSL "https://github.com/GZHoffie/skiver/releases/download/${SKIVER}/skiver" -o /usr/local/bin/skiver \
    && chmod +x /usr/local/bin/skiver \
    && curl -fsSL "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2}/minimap2-${MINIMAP2}_x64-linux.tar.bz2" \
        | tar -xj -C /opt && ln -s /opt/minimap2-${MINIMAP2}_x64-linux/minimap2 /usr/local/bin/minimap2 \
    && git clone --depth 1 -b "${MINIBWA}" https://github.com/lh3/minibwa /opt/minibwa && make -C /opt/minibwa \
    && ln -s /opt/minibwa/minibwa /usr/local/bin/minibwa \
    && curl -LsSf https://astral.sh/uv/0.10.0/install.sh | env UV_INSTALL_DIR=/usr/local/bin UV_NO_MODIFY_PATH=1 sh \
    && apt-get purge -y git make gcc libc6-dev && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*
COPY pyproject.toml uv.lock .python-version README.md /app/
COPY src /app/src
ENV UV_PYTHON_INSTALL_DIR=/opt/python UV_PROJECT_ENVIRONMENT=/opt/venv
RUN cd /app && uv sync --locked --no-dev && rm -rf /root/.cache /app/src/*
ENV PATH=/opt/venv/bin:$PATH
