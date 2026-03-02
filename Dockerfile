FROM rust:1-bookworm AS builder

WORKDIR /app

COPY fxsplit/Cargo.toml fxsplit/Cargo.lock ./
COPY fxsplit/src ./src

RUN cargo build --release --locked

FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/fxsplit /usr/local/bin/fxsplit

CMD ["fxsplit", "--help"]
