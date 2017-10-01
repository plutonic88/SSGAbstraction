nohup jre1.8.0_144/bin/java -Xmx16000m -classpath /home/users/abasak/exp/cplex.jar -Djava.library.path=/home/users/abasak/exp -jar 1.jar 5 5 15 5 1 > out1.txt
sleep 10
nohup jre1.8.0_144/bin/java -Xmx16000m -classpath /home/users/abasak/exp/cplex.jar -Djava.library.path=/home/users/abasak/exp -jar 2.jar 8 8 30 8 2 > out2.txt
sleep 10
nohup jre1.8.0_144/bin/java -Xmx16000m -classpath /home/users/abasak/exp/cplex.jar -Djava.library.path=/home/users/abasak/exp -jar 3.jar 10 10 40 10 2 > out3.txt
sleep 10
nohup jre1.8.0_144/bin/java -Xmx16000m -classpath /home/users/abasak/exp/cplex.jar -Djava.library.path=/home/users/abasak/exp -jar 4.jar 14 14 50 14 2 > out4.txt
sleep 10




