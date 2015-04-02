#!/home/vagrant/.rvm/rubies/ruby-2.2.0/bin/ruby

require 'rake/clean'

def method_name
  (/`(.*)'/.match(caller.first) ? $1 : 'Unresolved_Method_Name')
end

# This regular expresion greedy matches a floating point number
RE = /([0123456789\.]+)\s/

DASSL = %w{dassla.f ddassl.f}
LIBDN = File.join "..","lib"
LIBFN = "libdassl.a"
SRCDN = File.join "..","src"

# $ rake clean
CLEAN.include("*.o","*.exe","*-das.dat","*-das-n.dat","Makefile")
# $ rake clobber
CLOBBER.include(LIBDN)

desc "Assure the lib directory exists"
directory LIBDN

desc "Make the dassl library in lib"
task :make => LIBDN do |t|
  begin
    objs  = DASSL.collect{|fn| fn.ext('.o')}.join(' ')
    Dir.chdir(SRCDN) do
      DASSL.each{|fn| sh "gfortran -c -g -O2 #{fn}"}
      sh "ar rcvo #{LIBFN} #{objs}"
      mv LIBFN, LIBDN
      sh "ranlib #{File.join(LIBDN,LIBFN)}"
      sh "rm #{objs}"
    end
  rescue Exception => e
    puts "ERROR [#{method_name}]: #{e.to_s} (#{e.class})"
    puts "BACKTRACE: \n" + e.backtrace.join("\n")
  end
end

# Example command line
# 
# $ rake compile[chemakzo.f]
# $ ./chemakzo.exe
# $ rake compile[hires.f]
# $ ./hires.exe
# $ rake compile[medakzo.f]
# $ ./medakzo.exe
# $ rake compile[pollu.f]
# $ ./pollu.exe
# $ rake compile[taband.f]
# $ ./taband.exe
# $ rake compile[tafull.f]
# $ ./tafull.exe
# $ rake compile[tafulllc.f]
# $ ./tafulllc.exe
# $ rake compile[tba.f]
# $ ./tba.exe
desc "Compile a solver"
task :compile, [:example] do |t,args|
  args.with_defaults(:example => "tafulllc.f")
  puts "Args used: #{args}"
  makefile = File.read("RBMakefile").gsub(/##MYEXAMPLE##/,args[:example].ext(''))
  fn = "Makefile"
  File.write(fn, makefile)
  sh "make"
end

# Each example numerically integrated via a solver has a set
# of data files. The data files were created in a rake task.
#
# Each data file contains "run characteristics" as defined in
# 
# http://www.dm.uniba.it/~testset/report/testset.pdf
#
# We monitor and report the run characteristics
#
# Solver, RTOL, SCD, Total Steps, Accepted Steps, #F, #JAC, CPU
#
# Here SCD is the number of Significan Correct digits and we always
# set ATOL=RTOL*1.0E-2.
#
# Example command line
# 
# $ rake example[hires,"-32.0","-9.0"]
# $ rake example[chemakzo,"-32.0","-9.0"]
# $ rake example[medakzo,"-32.0","-9.0"]
# $ rake example[pollu,"-32.0","-9.0"]
# $ rake example[tba,"-12.0","-9.0"]
# $ rake example[transamp,"-27.7","-9.0"]
#
# This rake expects the examples to be compiled.
desc "Run Example"
task :run, [:ex, :minexp, :maxexp] do |t, args|
  # These arguments must be adjusted for a problem and solver to
  # cover an interval of RTOLs on which the solver can 
  # integrate the problem.
  args.with_defaults(:ex => "hires", :minexp => -32.0, :maxexp => -9.0)
  puts "Args used: #{args}"
  # Discretize the interval into npt equidistant mesh
  npt = 19
  # The stepsize on the mesh
  step = -(args[:maxexp].to_f-args[:minexp].to_f)/npt
  # Loop to run the problem 11 times and get statistical median
  (1..11).each do |nrun|
    fn = "#{args[:ex]}-#{"%02d" % nrun}.dat"
    sh "rm #{fn}" if File.exists?(fn)
    # Loop over the mesh
    (0..npt).each do |n|
      # Define the tolerances RTOL and ATOL
      atol = Math.exp(args[:maxexp].to_f + n*step)
      rtol = atol*100.0
      # A clever trick to input into running program
      # Append output after each run
      sh "{ echo #{rtol}; echo #{atol} ; } | ./#{args[:ex]}.exe >> #{fn}"
    end
  end
end

# Example command line
# 
# $ rake datamine[hires] > hires-mined.dat
# $ rake datamine[chemakzo] > chemakzo-mined.dat
# $ rake datamine[medakzo] > medakzo-mined.dat
# $ rake datamine[pollu] > polly-mined.dat
# $ rake datamine[tba] > tba-mined.dat
# $ rake datamine[transamp] > transamp-mined.dat
desc "Datamine the output files"
task :datamine, [:ex] do |t,args|
  tols = {}
  tol = nil
  arr = nil
  args.with_defaults(:ex => "hires")
  puts "Args used: #{args}"
  # Search in local directory for all files of the form [EXAMPLE-*.dat].
  #
  # Open each file and look for lines containing
  # 
  # RTol
  # scd of Y
  # number of integration steps
  # number of accepted steps
  # number of f evaluations
  # number of Jacobian evaluations
  # CPU-time used
  #
  # with exactly one possible match to a number. The number may
  # contain a decimal point.
  #
  # This code assumes multiple instances of the above ordered content!
  # A new instance of RTol triggers a new record. Multiple records (11)
  # with exactly the same RTol are expected to find a statistical
  # median.
  Dir.glob("#{args[:ex]}-*.dat").sort.each do |fn|
    File.readlines(fn).each do |line|
      if line.include? "RTol"
        tol = line.match("=(.*)$")[1].strip
        tols[tol] ||= []
        arr = [tol]
      end
      if line.include? "scd of Y"
        arr << line.split(/ /).uniq.collect{|e| e.strip}[-2]
      end
      if line.include? "Number of integration steps"
        arr << line.match(RE)[1].to_i
      end
      if line.include? "Number of accepted steps"
        arr << line.match(RE)[1].to_i
      end
      if line.include? "Number of f evaluations"
        arr << line.match(RE)[1].to_i
      end
      if line.include? "Number of Jacobian evaluations"
        arr << line.match(RE)[1].to_i
      end
      if line.include? "CPU-time used"
        arr << line.match(RE)[1].to_f
        tols[tol] << arr
        # Can print each record here as it is generated
        # puts arr.join(', ')
      end
    end
  end
  # The matches are ordered based on the key RTol.
  #
  # Visually inspect the output to make sure all run characteristics
  # are identical for each particular RTol!
  #
  # The SCD or CPU time is the only run characteristic which is
  # expected to vary. We extract the middle value or statistical
  # median.
  puts "All run characteristics for #{args[:ex]} sorted by RTol"
  puts "Visually inspect this output to make sure SCD only varies"
  arr = []
  tols.keys.sort_by{|k| k.to_f}.each do |key|
    scd = []
    tols[key].each do |a|
      # Print each record sorted by RTol
      puts a.join(', ')
      # Collect all SCD for a RTol value 
      scd << a[-1]
    end
    # Sort the CPU-times and select the middle
    # value of 11 expected values.
    v = scd.sort
    arr << (tols[key][-1] << v[6]).join(', ')
  end
  puts "The run characteristics for #{args[:ex]}"
  # Print the run characteristics
  puts arr
end

