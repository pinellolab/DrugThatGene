from DrugBuddyWeb import app
import subprocess as sb

sb.call('sudo /sbin/iptables -I INPUT 1 -i wlan0 -p tcp --dport 8080 -j ACCEPT',shell=True)
sb.call('sudo /sbin/iptables -I INPUT 1 -i wlan0 -p tcp --dport 22 -j ACCEPT',shell=True)
app.run(host='0.0.0.0', port=8080, debug=False)

