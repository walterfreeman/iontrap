foreach (@ARGV)
{
  print ("running: python parsepng.py $_\n");
  `python parsepng.py $_`;
}
