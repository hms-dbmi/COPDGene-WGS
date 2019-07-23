
gcloud dataproc clusters create linds \
--project=copd-genes \
--image-version=1.4-debian9 \
--master-machine-type=n1-standard-2 \
--worker-machine-type=n1-standard-8 \
--num-workers=2 --num-preemptible-workers=0 \
--network=ines-hail \
--metadata=JAR=gs://hail-common/builds/0.2/jars/hail-0.2-652d93ae3419bf45c69ea48f9a12cd96613e8119-Spark-2.4.0.jar,ZIP=gs://hail-common/builds/0.2/python/hail-0.2-652d93ae3419bf45c69ea48f9a12cd96613e8119.zip,MINICONDA_VERSION=4.4.10 \
--initialization-actions=gs://dataproc-initialization-actions/conda/bootstrap-conda.sh,gs://script-cluster/script.py
