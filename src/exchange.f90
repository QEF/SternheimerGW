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
!> This module provides the routines to evaluate the exchange self-energy.
!!
!! In real space the exchange energy is defined as
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(r, r') = -\sum_{nq} f_{nk-q} 
!!   \frac{\psi_{nk-q}^\ast(r) \psi_{nk-q}(r')}{|r - r'|}~.
!! \f}
!! Here, the \f$f_{nk}\f$ and \f$\psi_{nk}\f$ are the occupation and the
!! eigenfunction for band \f$n\f$ and wave vector \f$k\f$, respectively.
!! We are interested in the quantity in reciprocal space, which is related to
!! the real space self-energy in the following way
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(G, G') = \Omega \sum_{r,r'} e^{i (k + G) r}
!!   \Sigma_{k}^{\text{x}}(r, r') e^{-i (k + G') r'}~.
!! \f}
!! To evaluate this quantity, we conveniently evaluate the Coulomb potential
!! in reciprocal space
!! \f{equation}{
!!   V_k(G) = \sum_{r} \frac{e^{i (k + G) r}}{r} = \frac{4\pi}{|k + G|^2}~.
!! \f}
!! For practical calculations, we need to truncate the integration, which is
!! handled in the truncation_module.
!!
!! We remind of the fact that the wave function contains a Bloch phase and
!! a lattice periodic part
!! \f{equation}{
!!   \psi_{nk}(r) = e^{ikr} u_{nk}(r)
!! \f}
!! Replacing the Coulomb potential by the Fourier transform resolves the
!! convolution of \f$r\f$ and \f$r'\f$. The product of exponential and lattice
!! periodic part of the wave function corresponds to the Fourier transform of
!! the lattice periodic part. We obtain the following single \f$G''\f$ summation
!! \f{equation}{
!!   \Sigma_{k}^{\text{x}}(G, G') = -\sum_{nq} f_{nk-q} \sum_{G''}
!!   V_q(G'') u_{nk-q}^\ast(G - G'') u_{nk-q}(G'' - G')~.
!! \f}
MODULE exchange_module

  IMPLICIT NONE

CONTAINS

END MODULE exchange_module
