#!/usr/bin/env bash
#
# Written by Danilo de P. Kuritza <danilokuritza@gmail.com>
# Last updated on, Dec-18-2020

export LANG="en_US"
export LC_NUMERIC="en_US.UTF-8"

set -o errexit
set -o nounset

####################################################################
#    This script was created in order to calculate the biaxial     # 
#    strain (X and Y directions) in 2D materials. To work, copy    #
#    the POTCAR and POSCAR files of the desired material in the    #
#    same directory as ./strain.sh.                                #
#                                                                  #
#    How to run: nohup bash strain.sh > saida.out 2> job.err &     #
#                                                                  #
#    WARNING: If possible, install ssmtp via apt install before    #
#    runing the script!                                            #
####################################################################

###################
#    Variables    #
###################

### User setup ###
EMAIL="danilokuritza@gmail.com"        # User e-mail
MACHINE=`hostname`                     # Current machine
CDIR=`pwd`                             # Current directory
FLAG_MAIL="true"                       # true (send email) | false (don't send email)

### VASP setup ###
VASP_DIR=~/Tools/VASP/5.4.4/bin/
OUTPUT="job"                           # Output file name
NPROC="8"                              # Number of cores
NPAR="2"                               # Sets the NPAR number in INCAR file
KPAR="4"                               # Sets the KPAR number in INCAR file
KP="8"                                 # KPOINTS (Do not use large numbers if FHYB=true!) 
ENCUT="400"                            # Cutoff energy (eV)
FHYB="false"                           # Hybrid functionals (HSE06): true | false
VDWT="11"                              # vdW correction tag: 11 (DFT-D3 Zero) | 12 (DFT-D3 BJ) | 20 (TS) | 0 (OFF)
VDWSCS=".FALSE."                       # Self-consistent screening (use only if VDWT=20 !): .TRUE. | .FALSE.

### Strain setup ###
STEP="0.50"                            # Step used in the strain
NMIN="-1.00"                           # initial %
NMAX="1.00"                            # final %

################################
#    User Defined Functions    #
################################

incar_maker (){
if [ ${FHYB} = true ]; then
  LHF=".TRUE."
  PRE="F"
  AEX="0.25"
  HFS="0.2"
  ALG="All"
  ISY="-1"
else
  LHF=".FALSE."
  PRE="N"
  AEX="0"
  HFS="0.0"
  ALG="Normal"
  ISY="-1"
fi

cat > INCAR <<EOF 
PREC        = Accurate
LREAL       = .FALSE.
LZEROSTRESS = $1
ENCUT       = $2
ISMEAR      = 0
SIGMA       = 0.1
ALGO        = ${ALG}
ISIF        = 3
IBRION      = 1
ISYM        = ${ISY}
SYMPREC     = 1.0e-09
NSW         = 100
EDIFF       = 0.5e-07
EDIFG       = -1.0e-04
IVDW        = ${VDWT}
LVDWSCS     = ${VDWSCS}
LHFCALC     = ${LHF}
PRECFOCK    = ${PRE}
AEXX        = ${AEX}
HFSCREEN    = ${HFS}
LDIPOL      = .TRUE.
IDIPOL      = 3
NPAR        = ${NPAR}
KPAR        = ${KPAR}
EOF
}

kpoints_maker (){
cat > KPOINTS <<EOF
Automatic mesh
0
Gamma
$1  $1  1
0  0  0
EOF
}

### Inputs: $stnType $per $parx|$pary ###
poscar_modifier (){ 
  VAR=`echo "$2*(1+($1*0.01))" | bc -l`
  VARCUT=`printf %.16f "${VAR}"`
  awk -v VAR_IN=$2 -v VAR_OUT=${VARCUT} '{sub(VAR_IN,VAR_OUT)}1' POSCAR > POSCAR.temp && mv POSCAR.temp POSCAR
}

### Inputs: nothing!(global variables) ###
vasp_run (){
  ### See if VASP DIR is valid ###
  [ -d ${VASP_DIR} ] || error_exit "Unknown ERROR detected! Check the VASP directory (${VASP_DIR})."
  nohup mpirun -np "${NPROC}" "${VASP_DIR}"vasp_std > "${OUTPUT}".out 2> "${OUTPUT}".err || error_exit "Unknown ERROR detected! Check the VASP input and output files for more details."
  if cat ${OUTPUT}.out | grep -q "reached required accuracy"; then
    printf "Finished!\n"
  else
    error_exit "The calculation did not converge! Change the values for NSW, EDIFF and EDIFFG for better accuracy and try again."
  fi
}

strain_grabb (){
  parxd=`awk '{if(NR==3) print $1}' CONTCAR`
  paryd=`awk '{if(NR==4) print $2}' CONTCAR`
  epsx=`echo "(${parxd}-${parx0})/${parx0}" | bc -l`
  epsy=`echo "(${paryd}-${pary0})/${pary0}" | bc -l`
  strain=${epsx}
}

gap_grabb (){
  homo=`awk '/NELECT/ {print $3/2}' OUTCAR`
  lumo=`awk '/NELECT/ {print $3/2+1}' OUTCAR`
  nkpt=`awk '/NKPTS/ {print $4}' OUTCAR`
  e1=`grep "     ${homo}     " OUTCAR | head -${nkpt} | sort -n -k 2 | tail -1 | awk '{print $2}'`
  e2=`grep "     ${lumo}     " OUTCAR | head -${nkpt} | sort -n -k 2 | head -1 | awk '{print $2}'`
  if [ 1 -eq "$(echo "${e1} >= 0" | bc)" ]; then
    gap=`echo "sqrt(${e2}^2) - sqrt(${e1}^2)" | bc`
  elif [ 1 -eq "$(echo "${e2} <= 0" | bc)" ]; then
    gap=`echo "sqrt(${e1}^2) - sqrt(${e2}^2)" | bc`
  else
    gap=`echo "sqrt(${e1}^2) + sqrt(${e2}^2)" | bc`
  fi
}

energy_grabb (){
  energy=`grep 'without entropy=' OUTCAR | tail -1 | awk '{print $7}'`
}

data_printer (){
  printf "%5.2f %10.8f %10.8f %24.20f %12.8f %8.4f %8.4f %8.4f\n" $1 ${parxd} ${paryd} ${strain} ${energy} ${e1} ${e2} ${gap} >> ../$2.dat
}

run_function (){
  vasp_run
  strain_grabb
  energy_grabb
  gap_grabb
  data_printer ${per} "xy"
  (head -n 1 ../xy.dat && tail -n +2 ../xy.dat | sort -k 1 -n) > ../xy.temp
  mv ../xy.temp ../xy.dat
}

prep_relax (){
  mkdir ./rlx/ && cd ./rlx/ && cp ../{POSCAR,POTCAR} .
  incar_maker ".FALSE. .FALSE. .TRUE." ${ENCUT}
  kpoints_maker ${KP}
  vasp_run
  cd ../ && mv POSCAR POSCAR-old && cp ./rlx/CONTCAR ./POSCAR
}

error_exit(){
  echo "$1" 1>&2
  exit 1
}

######################
#    Main Program    #
######################

printf "
####################################################################
#    Written by Danilo de P. Kuritza <danilokuritza@gmail.com>     #
#    Last updated on, Dec-10-2020                                  #
#                                                                  #
#    This script was created in order to calculate the biaxial     # 
#    strain (X and Y directions) in 2D materials. To work, copy    #
#    the POTCAR and POSCAR files of the desired material in the    #
#    same directory as ./strain.sh.                                #
#                                                                  #
#    How to run: nohup bash strain.sh > saida.out 2> job.err &     #
#                                                                  #
#    WARNING: If possible, install ssmtp via apt install before    #
#    runing the script!                                            #
####################################################################\n\n"

START_TIME=$(date +%s)

### See if $file exists or not. If one of the files is not present, then exit status # 1 ###
printf "Checking POSCAR and POTCAR files: "
for file in POSCAR POTCAR; do
  [ ! -s ${file} ] && error_exit "WARNING: ${file} not found! You MUST place the POTCAR and POSCAR files in this same directory."
done
printf "Files are OK! Moving on...\n\n"

###  Initial relaxation ###
printf "Preparing structure: "
if [ -s POSCAR-old ] && [ -d ./rlx/ ] && [ -s ./rlx/CONTCAR ] && cmp -s ./POSCAR ./rlx/CONTCAR; then
  printf "Finished!\n\n"
else      
  if [ -d ./rlx/ ]; then
    rm -r ./rlx/
  fi
  if [ -s POSCAR-old ]; then
    mv POSCAR-old POSCAR
  fi
  prep_relax
  echo ""
fi

parx0=`awk '{if(NR==3) print $1}' POSCAR`
pary0=`awk '{if(NR==4) print $2}' POSCAR`

if [ -d "xy" ]; then
  flagRest="false"
else
  flagRest="true"
fi

if [ -d xy ]; then
  cd xy 
else 
  mkdir xy && cd xy
fi

### Replaces the string "m" with "-" ###
if [ "$(ls -A ./)" ]; then
  dirlist=(`ls -d ./*/`)
  for dir in ${dirlist[*]}; do
    if [[ ${dir} == *"m"* ]]; then
      mv -- "${dir}" "${dir//m/-}" 
    fi
  done
fi

if [ -s xy.dat ]; then
  printf "xy.dat already exists! Appending data...\n\n"
else
        printf "%5s %10s %10s %24s %12s %8s %8s %8s\n" % x y Strain Energy HOMO LUMO gap >> xy.dat
fi
for per in $(seq ${NMIN} ${STEP} ${NMAX}); do
  printf "Running ${per}%% strain: "
  if [ -d ./${per}/ ]; then
    cd ./${per}/ && cp ../../{POSCAR,POTCAR} .
  else
    mkdir ./${per}/ && cd ./${per}/ && cp ../../{POSCAR,POTCAR} .
  fi
  LZS=".TRUE. .TRUE. .TRUE."
  poscar_modifier ${per} ${parx0}
  poscar_modifier ${per} ${pary0}
  incar_maker "${LZS}" ${ENCUT}
  kpoints_maker ${KP}
  if [ "${flagRest}" = true ]; then
    run_function
  else
    if [ -f ${OUTPUT}.out ] && [ -s CONTCAR ]; then
      if cat ${OUTPUT}.out | grep -q "reached required accuracy"; then
        echo "Calculation already converged. Going to the next..."
      else
        rm CHG* EIGENVAL OUTCAR REPORT vasprun.xml XDATCAR DOSCAR IBZKPT ${OUTPUT}* OSZICAR PCDAT WAVECAR
        run_function
      fi
    else
      run_function
    fi
  fi
  outputList="POTCAR XDATCAR EIGENVAL WAVECAR IBZKPT OSZICAR vasprun.xml CHG* REPORT PCDAT DOSCAR"
  for output in ${outputList}; do
    if [ -f ${output} ]; then
      rm  ${output}
    fi
  done
  END_TIME=$(date +%s)
  printf "time: $((${END_TIME} - ${START_TIME})) s\n\n"
  cd ../
done

 ### Replaces the string "-" with "m" ###
 for dir in *"-"*; do
   mv -- "${dir}" "${dir//-/m}"
 done

### E-mail section (ssmtp) ###
if [ "${FLAG_MAIL}" = false ]; then
  printf "\nThis is the end my friend!\n"
else
  if printf "TO: ${EMAIL}\nSUBJECT: Strain Calculation\n\nThe strain calculation for your material is finally over! (MACHINE: ${MACHINE} - DIR: ${CDIR})" | ssmtp ${EMAIL}; then
    printf "email sent to: ${EMAIL}\n\n\nThis is the end my friend!\n"
  else
    error_exit "Cannot seend e-mail! See if ssmtp is properly setup."
  fi
fi

exit 0
