# ---------- Build Stage ----------
FROM rust:1-bookworm AS builder

WORKDIR /app

COPY fxsplit/Cargo.toml fxsplit/Cargo.lock ./
COPY fxsplit/src ./src

RUN cargo build --release --locked && \
    strip target/release/fxsplit

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/fxsplit /usr/local/bin/fxsplit

# Set up non-root user
RUN useradd -m -u 1000 fxsplituser && \
    chmod +x /usr/local/bin/fxsplit

USER fxsplituser
WORKDIR /data

RUN fxsplit --help

CMD ["fxsplit", "--help"]
