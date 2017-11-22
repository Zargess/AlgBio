
$ourtimes = New-Object System.Collections.Generic.List[System.Object]
$rapidtimes = New-Object System.Collections.Generic.List[System.Object]
$quicktimes = New-Object System.Collections.Generic.List[System.Object]
$files = Get-ChildItem "C:\Users\mathias\AlgBio\Project5\test" 
foreach($f in $files) {
	$outname1 =  $f.BaseName + "_ours.newick"
	$outname2 =  $f.BaseName + "_rapid.newick"
	$outname3 =  $f.BaseName + "_quick.newick"
	$sw = [Diagnostics.Stopwatch]::StartNew()
	python nj.py $f.FullName | Out-File $outname1 -encoding ASCII
	$sw.Stop()
	$ourtimes.Add($f.BaseName)
	$ourtimes.Add($sw.Elapsed)
	$sw = [Diagnostics.Stopwatch]::StartNew()
	.\rapidnj.exe -i pd $f.FullName | %{$_ -replace "'",""} | Out-File $outname2 -Encoding ASCII
	$sw.Stop()
	$rapidtimes.Add($f.BaseName)
	$rapidtimes.Add($sw.Elapsed)
	$sw = [Diagnostics.Stopwatch]::StartNew()
	.\quicktree.exe -in m $f.FullName | Out-File $outname3 -encoding ASCII
	$sw.Stop()
	$quicktimes.Add($f.BaseName)
	$quicktimes.Add($sw.Elapsed)
}
$ourtimes | Out-File "ourtimes.txt"
$rapidtimes | Out-File "rapidtimes.txt"
$quicktimes | Out-File "quicktimes.txt"