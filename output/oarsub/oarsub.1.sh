#!/bin/sh
#OAR -O ./stdout/oar.1.out
#OAR -E ./stdout/oar.1.err
#OAR --resource core=2,walltime=00:06:00
exec mpirun -machinefile "$OAR_NODE_FILE" \
  -- env MELISSA_JOB_ID=0 PATH=/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-melissa-core-develop-v4kaesy527uxmla3ikqqafauc2e6lmou/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-plotext-5.2.8-fpfvh57mejmva3or6ozivbdwjjimtxag/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-jsonschema-4.17.3-3lh53o7owyqgo6ptijquhroxeocgna2k/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-numpy-1.26.2-cggjy5fgiczh5fv5maeqowkutfjz7ofd/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/python-3.10.13-s2rx3ucgx7osywrxyxbrruxiuq66zg7a/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/xz-5.4.1-tgx4eljjy4cqn4efruveurnipahsfnct/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-melissa-core-develop-v4kaesy527uxmla3ikqqafauc2e6lmou/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-pybind11-2.11.0-g2gw7fqfm6zgf2y6xclc3x2hv7u2eokj/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libzmq-4.3.5-2uhcvgzinirkggmcse55sgeledo65xjb/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-plotext-5.2.8-fpfvh57mejmva3or6ozivbdwjjimtxag/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/openmpi-4.1.6-6ej7swthbz6qe6jfxwmtole7ufkormue/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/pmix-5.0.1-jho5uje2hwpulftrkqchehskmjjhfj2d/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libevent-2.1.12-uaqbnre2njocvcnod5533v3adlcfc4uz/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/numactl-2.0.14-2uwin2wnry5yam4rsyc4wh3u7itkqxut/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/hwloc-2.9.1-tpxdqo7c7suhezyo2npo7pvapurzovji/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-jsonschema-4.17.3-3lh53o7owyqgo6ptijquhroxeocgna2k/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/py-numpy-1.26.2-cggjy5fgiczh5fv5maeqowkutfjz7ofd/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/openblas-0.3.25-5bxupihx32yt5sxiukamli2hau6mrffj/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/python-3.10.13-s2rx3ucgx7osywrxyxbrruxiuq66zg7a/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/util-linux-uuid-2.38.1-53bvehimk7ttukdchwlw2pjhizmmyotz/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/sqlite-3.43.2-uykyxxlkk4ifpgxu6oltfwy4ubjmtb3v/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/openssl-3.1.3-pl2wnqd64ajfq7owu3pgauyohszbnye3/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/gettext-0.22.4-2f63iqendy7372kmvyjbbckmhlldbams/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libxml2-2.10.3-vf4hwqswuujjn6weri4qyfjq6l6gh6pt/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/xz-5.4.1-tgx4eljjy4cqn4efruveurnipahsfnct/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libiconv-1.17-cxccx7rvry5ptgiwyfyt33ldkgetvfqa/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/gdbm-1.23-74ua7wi4qmd3brgdnrvg3y7recqsulzx/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/readline-8.2-mlzdgeyt4obdumpltz3gtvmwjw4chixi/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/ncurses-6.4-vzpityyy2a54wajn6ipt6hhhhyzpykx2/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/expat-2.5.0-jnq4gpec5wiq3hbo2ubsketermkxpbmj/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/bzip2-1.0.8-3arpynddwetsswiniq2y37xbuyzurynp/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/melissa-api-develop-whvfbeandg277pfxc7sy6zwwmw633ibs/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/melissa-api-develop-whvfbeandg277pfxc7sy6zwwmw633ibs/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/openmpi-4.1.6-6ej7swthbz6qe6jfxwmtole7ufkormue/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/pmix-5.0.1-jho5uje2hwpulftrkqchehskmjjhfj2d/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libevent-2.1.12-uaqbnre2njocvcnod5533v3adlcfc4uz/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/openssl-3.1.3-pl2wnqd64ajfq7owu3pgauyohszbnye3/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/numactl-2.0.14-2uwin2wnry5yam4rsyc4wh3u7itkqxut/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/hwloc-2.9.1-tpxdqo7c7suhezyo2npo7pvapurzovji/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/ncurses-6.4-vzpityyy2a54wajn6ipt6hhhhyzpykx2/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libxml2-2.10.3-vf4hwqswuujjn6weri4qyfjq6l6gh6pt/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/xz-5.4.1-tgx4eljjy4cqn4efruveurnipahsfnct/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libiconv-1.17-cxccx7rvry5ptgiwyfyt33ldkgetvfqa/bin:/home/catzarakis/spack/opt/spack/linux-debian11-haswell/gcc-10.2.1/libzmq-4.3.5-2uhcvgzinirkggmcse55sgeledo65xjb/bin:/home/catzarakis/spack/bin:/home/catzarakis/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/grid5000/code/bin:/opt/puppetlabs/bin:/home/catzarakis/.dotnet:/home/catzarakis/.dotnet/tools ./client_scripts/client.0.sh
