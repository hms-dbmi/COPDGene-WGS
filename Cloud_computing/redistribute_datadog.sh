
# Amazon Web Services: all nodes
for WORKERIP in `sudo grep -i privateip /mnt/var/lib/info/*.txt | sort -u | cut -d "\"" -f 2`
do
   
   # Install data dog in worker nodes
   ssh hadoop@${WORKERIP} "DD_API_KEY=e5c74a4d3c394e204f9256402ae42880 bash -c \"\$(curl -L https://raw.githubusercontent.com/DataDog/datadog-agent/master/cmd/agent/install_script.sh)\""
done



# Google Dataproc: single node
gcloud compute ssh lindsay1-m: --zone us-east1-c
DD_API_KEY=fb3aa58156d97c2629705958bc954f2a bash -c "$(curl -L https://raw.githubusercontent.com/DataDog/datadog-agent/master/cmd/agent/install_script.sh)"
