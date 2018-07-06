!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
!> This module keeps track of the SternheimerGW version number.
!!
!! We use semantic versioning as much as is reasonably possible for a code
!! relying on the occationally rapid changes in the Quantum Espresso package.
!! A detailed description can be found here http://semver.org/ .
!!
!! In addition, we make use of a 
!! <a href="http://nvie.com/posts/a-successful-git-branching-model"> Gitflow</a>
!! <a href="https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow">
!! Workflow</a>.
!!
!! Here is an outline of the general principles for naming versions:
!! - Whenever there is an incompatible change to the API, update the major
!!   version number. In practice this means that if you make significant
!!   changes to the interface of a subroutine or the input file.
!!   Significant are all those changes that would break an old file unless
!!   it is modified to the new behavior.
!! - Whenever a new feature is introduced and whenever QE is updated, update
!!   the minor revision number. Every time you do a minor release update the 
!!   compatible_QE_revision field to the correct revision, so that one knows
!!   with which revision of QE the code should compile.
!! - Use the patch sublevel only if you must fix the master between releases.
!!   Typically the QE release cycle should be so fast that this won't be
!!   necessary, but if a very severe bug occurs fix it and release a new
!!   version increasing the third digit.
!! - Use tags to identify releases and add the compatible QE release. Example:
!!   If version 1.15 of SternheimerGW is released with QE version 6.2 use the tag
!!   v1.15_QE6.2. Do only use this form of the tag if the version corresponds
!!   excactly to the unmodified version of QE; otherwise use simply the
!!   v1.15 tag. Note that we use '_' instead of '+', because Github replaces
!!   '+' with '-' in the filename. Then the download script of QE does not work.
!! - For versions between commits the git-describe functionality is used to
!!   determine the version number.
!!
!! Here is an outline of the Gitflow Workflow
!! - The master branch of the project should be in production ready state all
!!   the time. You want to update the master typically with every release of
!!   QE so that SternheimerGW remains compatible and occassionally in between if an
!!   important feature is added and the next QE release is far away.
!! - To create a new release fork a release branch from the develop branch and
!!   do all the necessary adjustments there. This includes the update of the
!!   version number, bugfixes if necessary, and any last minute adjustment if
!!   the QE API changes. When the release branch is ready merge it into the
!!   master using --no-ff to create a commit. Also merge it back into develop
!!   so that the version number is updated there, too. There **must** not be any
!!   new feature added to the release branch after it forked. Check that the
!!   branch passes the buildbot before merging it into the master!
!! - The develop branch is the working branch in which new features are added
!!   for a future release. It is supposed to be in a working state at all time,
!!   but occassionally it might be broken. For this the branch is constantly
!!   monitored by the buildbot. On this branch you should see the following
!!   commits: hotfixes to master (merged into develop); bugfixes due to API
!!   changes in QE; bugfixes to the develop branch detected by the devloper or
!!   the buildbot; and importantly new feature commits (see below). Keep in mind
!!   to update the QE compatibility when you commit to develop.
!! - The main body of your work should take place in feature branches. In these
!!   you can play around and try out new stuff. If you develop a complex feature
!!   try to see if you can chop it into smaller logical units, so that the
!!   lifetime of a feature branch doesn't get to long. When you are confident
!!   that the feature is stable you merge it into the develop branch so that it
!!   get's roled out into the next release. If you find someone to code review
!!   your changes before merging all the better. When you merge the feature
!!   into develop, remember to update the release notes, so that your changes
!!   get documented in the next release.
!! - On the public repository only the master and the develop branch are stored.
!!   If you want to share a feature branch, use a local repository instead.
!!
!----------------------------------------------------------------------------
MODULE gw_version
  !
  IMPLICIT NONE
  !
  !> Version number set by hand.
  CHARACTER(*), PARAMETER :: gw_version_number = '0.15'
  !> Version number extracted from git describe.
  CHARACTER(*), PARAMETER :: gw_git_describe = 'unknown'
  !> Git commit of QE compatible with the commit of SternheimerGW
  CHARACTER(*), PARAMETER :: compatible_QE_commit = 'QE 6.3 7357cdb'
  !
END MODULE gw_version
