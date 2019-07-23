
# https://aws.amazon.com/premiumsupport/knowledge-center/emr-timeout-connection-wait/
# In linux connect to master node
# ssh -i ~/.ssh/hail_avillachlab_73_key.pem hadoop@ec2-3-92-22-64.compute-1.amazonaws.com

# when in master node:
sudo vi /usr/share/aws/emr/emrfs/conf/emrfs-site.xml

#change it to 1000 for pruning:
<configuration>
  <property>
  <name>fs.s3.maxConnections</name>
  <value>1000</value>
  </property>
</configuration>


#sudo cp emrfs-site.xml /usr/share/aws/emr/emrfs/conf/emrfs-site.xml

for WORKERIP in `sudo grep -i privateip /mnt/var/lib/info/*.txt | sort -u | cut -d "\"" -f 2`
  do
     # Distribute keys to workers for account hadoop
     scp -o "StrictHostKeyChecking no" /usr/share/aws/emr/emrfs/conf/emrfs-site.xml ${WORKERIP}:/tmp/emrfs-site.xml
     ssh hadoop@${WORKERIP} "sudo cp /tmp/emrfs-site.xml /usr/share/aws/emr/emrfs/conf/emrfs-site.xml"
     sudo stop hadoop-yarn-resourcemanager; sudo start hadoop-yarn-resourcemanager
done

#restart hadoop
sudo stop hadoop-yarn-resourcemanager; sleep 1; sudo start hadoop-yarn-resourcemanager
