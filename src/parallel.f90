!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
!> This module provides the subroutines to parallelize a task over a communicator.
MODULE parallel_module

  IMPLICIT NONE

CONTAINS

  !> Parallelize a given number of tasks over given communicator.
  !!
  !! Distribution of tasks over processes in communicator according to the
  !! following strategy
  !!
  SUBROUTINE parallel_task(comm, num_task_total, first_task, last_task, num_task)

    USE mp_global, ONLY: mp_rank, mp_size

    !> the communicator over which the tasks are distributed
    INTEGER, INTENT(IN)  :: comm

    !> total number of tasks to assign
    INTEGER, INTENT(IN)  :: num_task_total

    !> first task assigned to this process
    INTEGER, INTENT(OUT) :: first_task

    !> last task assigned to this process
    INTEGER, INTENT(OUT) :: last_task

    !> the number of tasks assigned to each process
    !! @note you can use this for the mp_gather wrapper in this module
    INTEGER, INTENT(OUT), ALLOCATABLE :: num_task(:)

    !> number of processes in communicator
    INTEGER num_proc

    !> rank of this process in communicator
    INTEGER my_rank

    !> minimal number of task per process
    INTEGER num_task_min

    !> remainder after assigning an equal amount of processes
    INTEGER num_remain

    !> at this point we add one additional task
    INTEGER last_proc

    ! determine rank and size
    ! note: add 1 to rank, because Fortran indices start counting at 1
    my_rank = mp_rank(comm) + 1
    num_proc = mp_size(comm)

    ! allocate array for assigned tasks
    ALLOCATE(num_task(num_proc))

    !! 1. Distribute num_task / size(comm) on every process.
    num_task_min = num_task_total / num_proc

    !! 2. Determine the remainder of num_remain of unassigned tasks.
    num_remain = MOD(num_task_total, num_proc)

    !! 3. Assign the last processes in the list an extra task.
    last_proc = num_proc - num_remain
    num_task(:last_proc) = num_task_min
    num_task(last_proc+1:) = num_task_min + 1

    !! 4. Determine first and last task for current process.
    last_task = SUM(num_task(:my_rank))
    first_task = last_task - num_task(my_rank) + 1

END SUBROUTINE parallel_task

END MODULE parallel_module
