Last login: Thu Feb 11 10:27:07 on ttys003
(base) sarehchimeh@Sarehs-MacBook-Air ~ % cd 
(base) sarehchimeh@Sarehs-MacBook-Air ~ % cd PycharmProjects 
(base) sarehchimeh@Sarehs-MacBook-Air PycharmProjects % ls
Surfaces	git		pythonProject	pythonProject1
(base) sarehchimeh@Sarehs-MacBook-Air PycharmProjects % cd Surfaces 
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % ls
Downsizing		Viruses			practice.py
README.md		overlapping_genes	retrieve_sequences
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote -v 
origin	https://github.com/Sarehchimeh/Surfaces.git (fetch)
origin	https://github.com/Sarehchimeh/Surfaces.git (push)
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote set-url poonlab https://github.com/PoonLab/Surfaces.git
fatal: No such remote 'poonlab'
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote add poonlab https://github.com/PoonLab/Surfaces.git
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git remote -v
origin	https://github.com/Sarehchimeh/Surfaces.git (fetch)
origin	https://github.com/Sarehchimeh/Surfaces.git (push)
poonlab	https://github.com/PoonLab/Surfaces.git (fetch)
poonlab	https://github.com/PoonLab/Surfaces.git (push)
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % git push poonlab master
error: src refspec master does not match any
error: failed to push some refs to 'https://github.com/PoonLab/Surfaces.git'
(base) sarehchimeh@Sarehs-MacBook-Air Surfaces % ssh sareh@129.100.26.213 
sareh@129.100.26.213's password: 
Permission denied, please try again.
sareh@129.100.26.213's password: 
Welcome to Ubuntu 18.04.5 LTS (GNU/Linux 4.15.0-132-generic x86_64)

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/advantage

 * Canonical Livepatch is available for installation.
   - Reduce system reboots and improve kernel security. Activate at:
     https://ubuntu.com/livepatch

58 packages can be updated.
3 updates are security updates.

New release '20.04.2 LTS' available.
Run 'do-release-upgrade' to upgrade to it.

*** System restart required ***
Last login: Wed Feb 10 20:36:38 2021 from 99.251.234.133
sareh@Erasmas:~$ cd script/
sareh@Erasmas:~/script$ ls
fasttree_iterator.sh       prune_CDS.py           prunetree_iterator.sh
gotoh2                     pruned_accn.py         prunetree.py
performance_prune_CDS.txt  prunetree_iterator.py  test.sh
sareh@Erasmas:~/script$ nano prune_CDS.py 

  GNU nano 2.9.3                     prune_CDS.py                               

                    for i in x:
                        if i in accessions:
                            pruned_record_count +=1
                            SeqIO.write(record,out_file,'fasta')
        print('record count: '+ str(record_count))
        print('pruned record count: '+ str(pruned_record_count))
                                    # write the next line into a file (record s$

if __name__ == "__main__":
    main()
