
## initialize a container named blendit in background
```
docker run -t -d --rm 
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro 
           -u $(id -u $USER):$(id -g $USER) 
           -v $PWD:/mnt 
           --name blendit shengwei/blendit:latest
```

## run blendit 
```
# profile kmer 
cmd1="blendit profile kmer [OPTIONS] ASSEMBLY"
echo $cmd1
docker exec dastool /bin/bash -c "$cmd1"

## kill container
docker container kill blendit

