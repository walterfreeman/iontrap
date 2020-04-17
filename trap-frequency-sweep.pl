#<this> <ac frequency in GHz> <az_min> <az_max> <az_step> <qz_min> <qz_max> <qz_step>\n");
#    exit(1);

for ($w = 10; $w>0.1; $w *= 0.9)
{
  $omegalabel = sprintf("%.2fGHz",$w);
  system("./trap $w -10 10 0.05 0 10 0.05 | grep RES > iontrap-results-$omegalabel-1e2cycles");
}
