# 1. Configure cluster

## 1.1. Get required software

```

sudo apt update -y
sudo apt install -y openjdk-8-jre-headless subversion

cd /shared/ ; \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh ; \
bash Miniconda3-latest-Linux-x86_64.sh ; \
source ~/.bashrc ;\
rm Miniconda3-latest-Linux-x86_64.sh ; \
conda install -c conda-forge mamba

svn checkout --username=yassinesouilmi svn+ssh://yassinesouilmi@svn.code.sf.net/p/human-te-dynamics/svn/ human-te-dynamics-svn ; \
cd /shared/human-te-dynamics-svn/ ; \
mamba env update -n base --file environment.yaml 

```

## 1.2. Fix slurm

1) First log in to the compute node from the head node and grab the "RealMemory" config in the node by running /opt/slurm/sbin/slurmd -C (as you did already)
2) Log out of the compute node back to the head node
3) Edit /opt/slurm/etc/slurm.conf and add NodeName=DEFAULT RealMemory=[RealMemory number] BEFORE include slurm_parallelcluster_nodes.conf
4) Log back into the compute node 
5) Run sudo service slurmd restart
5) Log back out of the compute node 
7) Run sudo service slurmctl restart
8)Check your work with scontrol show nodes

NodeName=ip-10-255-1-17 CPUs=32 Boards=1 SocketsPerBoard=1 CoresPerSocket=16 ThreadsPerCore=2 RealMemory=127459
UpTime=0-00:16:01