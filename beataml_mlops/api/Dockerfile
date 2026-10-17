# Build context is this directory. Run `python3 prepare_artifacts.py` first so
# ./artifacts/ holds read-only copies of the frozen model files to COPY in.
FROM python:3.11-slim

WORKDIR /app

# workflow.py's setdefault() leaves this alone if we set it first; /private/tmp
# (its own default, inherited from the original macOS development machine) may
# not be writable here.
ENV MPLCONFIGDIR=/tmp/mpl
RUN mkdir -p /tmp/mpl

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY artifacts ./artifacts
COPY main.py .

RUN useradd --create-home --shell /bin/false apiuser \
    && chown -R apiuser:apiuser /app /tmp/mpl
USER apiuser

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health', timeout=3)" || exit 1

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
