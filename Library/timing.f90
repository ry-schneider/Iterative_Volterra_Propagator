!> Timing
!! Provides code for timing code segements and printing out the times.
!! This only enables timing a single section at a time
!!
!! Author: Henry J Schmale

module timing
  implicit none
  private
  real        :: cpu_ti, cpu_tf
  integer     :: sys_ti, sys_tf, clockrate

  public begin_timing, stop_timing, print_timing
contains

  !> Starts the timer
  subroutine begin_timing
    call cpu_time(cpu_ti)
    call system_clock(sys_ti, clockrate)
  end subroutine begin_timing

  !> Stops the timer
  subroutine stop_timing
    call cpu_time(cpu_tf)
    call system_clock(sys_tf, clockrate)
  end subroutine stop_timing

  !> Get the elapsed time in seconds for CPU time
  function getCpuTimeElapsed() result(t_elaps)
    real(8) :: t_elaps
    t_elaps = cpu_tf - cpu_ti
  end function getCpuTimeElapsed 

  !> Get the wall clock time elapsed in seconds
  function getSystemTimeElapsed() result(t_elaps)
    real(8) :: t_elaps
    t_elaps = real(sys_tf - sys_ti) / clockrate
  end function getSystemTimeElapsed 

  !> Prints the results of the timing to standard output
  subroutine print_timing
    print*, "System Time Elapsed: ", getSystemTimeElapsed()
    print*, "Cpu Time Elapsed: ", getCpuTimeElapsed()
  end subroutine print_timing
end module timing

