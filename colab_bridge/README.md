# Colab GPU Bridge

This bridge lets the local autoresearch loop keep control while a manually
started Colab notebook acts as a GPU worker for DeepSurv jobs.

The shared queue directory must be visible to both sides. In practice, use a
Google Drive folder:

- Colab side: `/content/drive/MyDrive/autosota_colab_queue`
- Local side: the matching Google Drive Desktop path, passed with `--queue-root`

Do not put `affyfRMATest.csv` in the Colab data directory or queue. The worker
copies only:

- `affyfRMATrain.csv`
- `affyfRMAValidation.csv`
- `LOOCV_Genes2.csv`

## Local Side

Submit a DeepSurv GPU job after committing the candidate:

```bash
python colab_bridge/local_colab_queue.py submit \
  --queue-root "/path/to/Google Drive/autosota_colab_queue" \
  --job-id deepsurv_iter_001 \
  --require-pushed \
  --env DEEPSURV_ARENA_N_TRIALS=20 \
  --env DEEPSURV_ARENA_BOOTSTRAPS=2 \
  --env DEEPSURV_ARENA_EPOCHS=200
```

Watch and collect:

```bash
python colab_bridge/local_colab_queue.py watch --queue-root "/path/to/Google Drive/autosota_colab_queue" --job-id deepsurv_iter_001
python colab_bridge/local_colab_queue.py collect --queue-root "/path/to/Google Drive/autosota_colab_queue" --job-id deepsurv_iter_001
```

`collect` extracts the returned arena run artifacts into the local repo and
prints the result row.

## Colab Side

In Colab, mount Drive, verify GPU, then run:

```python
from google.colab import drive
drive.mount("/content/drive")

!python /content/drive/MyDrive/path/to/colab_gpu_worker.py \
  --queue-root /content/drive/MyDrive/autosota_colab_queue \
  --data-dir /content/drive/MyDrive/nsclc_train_valid_only \
  --repo-url https://github.com/osun24/autosota-nsclc-adj-chemo.git
```

The worker leases one queued job at a time. While a job runs, it updates a
heartbeat in `running/<job_id>.json`. If Colab disconnects, the heartbeat stops.
The next worker startup requeues stale running jobs after their lease expires.

## Disconnect Handling

- Jobs move `queued -> running -> done` or `failed`.
- Running jobs include `last_heartbeat` and `lease_expires_at`.
- A worker startup calls stale recovery before pickup.
- Local `requeue-stale` can also recover stale running jobs.
- Command failures go to `failed/` and are not retried automatically.
- Disconnects are retried by requeueing, up to `max_attempts`.
