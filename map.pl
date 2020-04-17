# fork some threads

$nthreads=64;
$created=0;
$cycles=100;

$wmax = 20;
$wmin = 0.1;
$wmult = 0.9;

$pid=1;
$threadid=0;
while ($created < $nthreads)
{
        if ($pid) {
                $created++;
                $pid=fork();
        }

        if ($pid)
        {
                $children[$created] = $pid;
        }
        else
        {
                $threadid=$created;
                last;
        }
}

#print "I am process ID $$ and my identity is $threadid\n";

for ($w = $wmax; $w > $wmin; $w*=$wmult ** $nthreads)
{
  $wtrue = $w * $wmult ** ($threadid);
  $wstring=sprintf("%05.3fGHz",$wtrue);
      $Vactrue = $Vac + $threadid * $vacstep;
     system("./trap $wtrue -10 10 0.05 0 10 0.05 1e-11 $cycles | grep RES > results/results-freq$wstring-cycles$cycles");
}

