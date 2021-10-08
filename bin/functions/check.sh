#!/bin/bash

check_path(){

    if [ ! -d "$HOME/bin" ]; then
        mkdir -p $HOME/bin
        cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
    export PATH="\$HOME/bin:\$PATH"
fi

EOF
    else
        perl -e 'exit(!(grep(m{^$ENV{HOME}/bin$},split(":", $ENV{PATH}))) > 0)'
        if [ $? != 0 ]; then
            cat <<EOF >> $HOME/.bashrc

if [ -d "\$HOME/bin" ] ; then
    export PATH="\$HOME/bin:\$PATH"
fi

EOF
        fi
    fi

. ~/.bashrc
}

check_connection(){
    echo -e "\nNetwork Status\n" &>> $log_file
    adress=$(ip r | grep default | cut -d ' ' -f 3)
    up=$(ping -q -w 1 -c 1 "$adress" > /dev/null && echo ok || echo error)
    if [ "$up" == "ok" ]; then
        echo -e "\nSESSION_TYPE = local\n" &>> $log_file
    elif [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ]; then
        echo -e "\nSESSION_TYPE = remote/ssh\n" &>> $log_file
    else
        echo -e "Network Error, please check your internet connection\n" &>> $log_file
        exit 1
    fi
}

check_dependencies(){
    counter=0
    printf '\n%s\t%20s\n' "DEPENDENCY" "STATUS"
    printf '%s\t%20s\n'   "----------" "------"

    for command in "$@"
    do
        length_command=$(echo $command | wc -m)
        distance_table=$((30 - $length_command))
        distance_expression=$(echo "%${distance_table}s")
        printf '%s' $command
        if ! [ -x "$(which $command 2> /dev/null)" ]; then
            printf $distance_expression
            printf "NOT INSTALLED\n"
            let counter++
        else
            printf $distance_expression
            printf "INSTALLED\n"
        fi
    done

    if [ $counter -gt 0 ]; then
        printf "ERROR: $counter missing dependencies, aborting execution\n" &>> $log_file
        exit 1
    fi
}

check_dockers(){
    counter=0
    printf '\n%s\t%28s\n' "DOCKER" "STATUS"
    printf '%s\t%28s\n' "------" "------"

    for docker in "$@"
    do
        length_docker=$(echo $docker | wc -m)
        distance_table=$((30 - $length_docker))
        distance_expression=$(echo "%${distance_table}s")
        printf '%s' $docker

        if [[ "$(docker images -q $docker 2> /dev/null)" == "" ]]; then
            printf $distance_expression
            printf "NOT INSTALLED\n"
            let counter++
        else
            printf $distance_expression
            printf "INSTALLED\n"
        fi
    done

    if [ $counter -gt 0 ]; then
        printf "ERROR: $counter missing dockers, aborting execution\n" &>> $log_file
        exit 1
    fi
}

check_databases(){

    dbNCBI="/mnt/disk1/bin/16S/NCBI.gz"
    dbRDP="/mnt/disk1/bin/16S/RDP.gz"
    dbSILVA="/mnt/disk1/bin/16S/SILVA.gz"
    dbKraken="/mnt/disk2/bin/Kraken2/yggdrasil"
    #Old version of Kmerfinder database
    # dbKmerfinder_prefix="/mnt/disk1/bin/kmerfinder_DB/bacteria.organisms.ATGAC"
    # dbKmerfinder_list=(.desc.p .p .len.p .ulen.p)
    dbKmerfinder_prefix="/mnt/disk1/bin/KmerFinder_DB/bacteria.ATG."
    dbKmerfinder_list=(seq.b name length.b comp.b)
    dbKmerfinder="${dbKmerfinder_list[@]/#/$dbKmerfinder_prefix}"
    databases="$dbNCBI $dbRDP $dbSILVA $dbKraken $dbKmerfinder"
    plasmid_db="/mnt/disk1/bin/plasmidid_db/plasmid.complete.nr100.fna"

    counter=0
    printf '\n%s\t%60s\n' "DATABASE" "STATUS"
    printf '%s\t%60s\n'   "--------" "------"

    for db in ${databases}
    do
        length_db=$(echo ${db} | wc -m)
        distance_table=$((70 - $length_db))
        distance_expression=$(echo "%${distance_table}s")
        printf '%s' $db

        if [ ! -e "$db" ] && [ ! -f  "$db" ]; then
            printf $distance_expression
            printf "NOT INSTALLED\n"
            let counter++
        elif [ ! -e "$db" ] && [ ! -d "$db" ]; then
            printf $distance_expression
            printf "NOT INSTALLED\n"
            let counter++
        else
            printf $distance_expression
            printf "INSTALLED\n"
        fi
    done

    if [ $counter -gt 0 ]; then
        printf "ERROR: $counter missing DB, aborting execution\n" &>> $log_file
        echo
        printf "INSTALL DATABASES: cd ottotype/ && bash databases.sh\n" &>> $log_file
        exit 1
    fi
}

checklist(){
    check_path &>> $log_file || error ${LINENO} $(basename $0)
    check_connection &>> $log_file || error ${LINENO} $(basename $0)
    check_dependencies salmonella.py minimap2 translate.py trimmomatic \
            sga Rscript idba_ud500 idba_ud spades.py bwa rename \
            samtools mlst kraken2 bioawk cd-hit &>> $log_file || error ${LINENO} $(basename $0)
    check_dockers seqsero srst2 ariba kraken2 kmerfinder plasmidid &>> $log_file || error ${LINENO} $(basename $0)
    check_databases &>> $log_file || error ${LINENO} $(basename $0)
}
