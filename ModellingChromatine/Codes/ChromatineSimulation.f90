program ChromatinSimulation
  implicit none

  ! Parameters for simulation
  integer, parameter :: chromatine_size = 60
  integer, parameter :: polymerase_count = 0
  integer, parameter :: simulation_steps = 5000
  integer, parameter :: adding_position = 25
  integer, parameter :: end_of_replication_position = chromatine_size - 25

  ! Simulation-specific parameters
  real(8), parameter :: histone_modification_percentage = 0.5d0
  real(8), parameter :: recruitment_probability = 1.0d0
  real(8), parameter :: alpha = 3.0d0 / 4.0d0
  real(8), parameter :: change_probability = alpha
  real(8), parameter :: regeneration_probability = 0.3d0
  real(8), parameter :: adding_polymerase_probability = 0.3d0
  real(8), parameter :: noisy_transition_probability = 1.0d0 - alpha
  integer, parameter :: vicinity_size = 5
  real(8), parameter :: F = alpha / (1.0d0 - alpha)

  ! Linear function parameters
  real(8), parameter :: slope = 1.0d-5
  real(8), parameter :: intercept = 0.0d0

  ! Polymerase movement probabilities
  real(8), parameter :: left_movement_probability = 0.5d0
  real(8), parameter :: right_movement_probability = 0.5d0

  ! Set seed for reproducibility
  call srand(42)

  type :: Chromatine
    character(1), dimension(:), allocatable :: histones
  contains
    procedure :: initialize => Chromatine_Initialize
    procedure :: noisy_transition => Chromatine_NoisyTransition
    procedure :: add_polymerases => Chromatine_AddPolymerases
    procedure :: adding_poly_proba => Chromatine_AddingPolyProba
    procedure :: change_next_histones => Chromatine_ChangeNextHistones
  end type Chromatine

  type :: Polymerase
    type(Chromatine) :: chromatine
    integer :: position
    real(8) :: temperature
  contains
    procedure :: move => Polymerase_Move
    procedure :: change_histones => Polymerase_ChangeHistones
    procedure :: delete => Polymerase_Delete
  end type Polymerase

  ! Global variables
  type(Chromatine) :: chromatine
  type(Polymerase), dimension(:), allocatable :: polymerases
  integer, dimension(:), allocatable :: existing_polymerase_positions
  integer :: polymerase_count_over_time(0:simulation_steps)
  integer :: active_histone_count_over_time(0:simulation_steps)
  integer :: acetylated_histone_count_over_time(0:simulation_steps)
  integer :: methylated_histone_count_over_time(0:simulation_steps)
  integer :: unmodified_histone_count_over_time(0:simulation_steps)
  integer :: noisy_changes_count_over_time(0:simulation_steps)
  integer :: enzyme_changes_count_over_time(0:simulation_steps)

  ! Open a file for data output
  open(unit=10, file='simulation_data.txt', status='replace')

  ! Initialize chromatine and polymerases with a specified temperature
  call chromatine%initialize(chromatine_size)
  call initialize_polymerases(existing_polymerase_positions, polymerase_count, adding_position)

  ! Perform simulation steps
  do frame = 0, simulation_steps
    ! Move and change histones for each polymerase
    do i = 1, size(polymerases)
      call polymerases(i)%move()
      call polymerases(i)%change_histones(change_probability, noisy_changes)
    end do

    ! Add new polymerases based on a probability
    if (rand() < adding_polymerase_probability) then
      call chromatine%add_polymerases(1, existing_polymerase_positions, adding_position)
    end if

    ! Delete polymerases that reach the end of replication
    do i = size(polymerases), 1, -1
      if (polymerases(i)%position >= end_of_replication_position) then
        call polymerases(i)%delete(existing_polymerase_positions)
        enzyme_changes = enzyme_changes + 1
      end if
    end do

    ! Update counts over time
    polymerase_count_over_time(frame) = size(existing_polymerase_positions) - 1
    active_histone_count_over_time(frame) = count(chromatine%histones == 'M' .or. chromatine%histones == 'A')
    acetylated_histone_count_over_time(frame) = count(chromatine%histones == 'A')
    methylated_histone_count_over_time(frame) = count(chromatine%histones == 'M')
    unmodified_histone_count_over_time(frame) = count(chromatine%histones == 'U')
    noisy_changes_count_over_time(frame) = noisy_changes
    enzyme_changes_count_over_time(frame) = enzyme_changes

    ! Write data to the file
    write(10, *) frame, polymerase_count_over_time(frame), active_histone_count_over_time(frame), &
                  acetylated_histone_count_over_time(frame), methylated_histone_count_over_time(frame), &
                  unmodified_histone_count_over_time(frame), noisy_changes_count_over_time(frame), &
                  enzyme_changes_count_over_time(frame)
  end do

  ! Close the data file
  close(10)

  contains

    subroutine Chromatine_Initialize(this, histones_count)
        class(Chromatine), intent(inout) :: this
        integer, intent(in) :: histones_count

        ! Initialize chromatine with histones
        allocate(this%histones(histones_count))
        this%histones = 'M'

        ! Randomly set approximately histone_modification_percentage of histones to 'U' (unmodified)
        this%histones(1:int(histone_modification_percentage * histones_count)) = 'U'
    end subroutine Chromatine_Initialize

    subroutine Chromatine_NoisyTransition(this, position, noisy_transition_probability, noisy_changes)
        class(Chromatine), intent(inout) :: this
        integer, intent(in) :: position
        real(8), intent(in) :: noisy_transition_probability
        integer, intent(inout) :: noisy_changes

        ! Implement the noisy transition logic
        if (rand() < noisy_transition_probability / 3.0d0) then
        if (this%histones(position) == 'A') then
            this%histones(position) = 'U'
            noisy_changes = noisy_changes + 1
        elseif (this%histones(position) == 'M') then
            this%histones(position) = 'U'
            noisy_changes = noisy_changes + 1
        elseif (this%histones(position) == 'U') then
            if (rand() < 0.5d0) then
            this%histones(position) = 'A'
            noisy_changes = noisy_changes + 1
            else
            this%histones(position) = 'M'
            noisy_changes = noisy_changes + 1
            end if
        end if
        end if
    end subroutine Chromatine_NoisyTransition

    subroutine Chromatine_AddPolymerases(this, count, existing_polymerase_positions, adding_position)
        class(Chromatine), intent(inout) :: this
        integer, intent(in) :: count
        integer, dimension(:), intent(inout) :: existing_polymerase_positions
        integer, intent(in) :: adding_position

        ! Add a specified number of new polymerases at non-overlapping positions
        integer :: new_position
        do i = 1, count
        new_position = adding_position
        do while (new_position .in. existing_polymerase_positions)
            new_position = new_position + 1
        end do
        existing_polymerase_positions = [existing_polymerase_positions, new_position]
        end do
    end subroutine Chromatine_AddPolymerases

    real(8) function Chromatine_AddingPolyProba(this, adding_position)
        class(Chromatine), intent(in) :: this
        integer, intent(in) :: adding_position

        ! Linear function of the local density of histones
        ! Let's calculate local density as the count of active histones in the vicinity of the polymerase
        integer :: start_index, end_index, local_density

        start_index = max(1, adding_position - vicinity_size)
        end_index = min(size(this%histones), adding_position + vicinity_size)
        local_density = count(this%histones(start_index:end_index) == 'M' .or. this%histones(start_index:end_index) == 'A')

        ! Calculate the probability of adding a new polymerase
        Chromatine_AddingPolyProba = slope * real(local_density) + intercept
    end function Chromatine_AddingPolyProba

    subroutine Chromatine_ChangeNextHistones(this, position, p_recruitment, p_change, enzyme_changes, nth_neighbor)
        class(Chromatine), intent(inout) :: this
        integer, intent(in) :: position, nth_neighbor
        real(8), intent(in) :: p_recruitment, p_change
        integer, intent(inout) :: enzyme_changes

        ! Simulate the influence of neighbors on the next histones
        if (1 <= position .and. position < size(this%histones) - 1) then
        character(1) :: current_histone
        integer :: nth_position

        current_histone = this%histones(position)

        ! Calculate the influence of vicinity on the recruitment probability
        real(8) :: adjusted_p_recruitment
        adjusted_p_recruitment = p_recruitment
        if (adjusted_p_recruitment > 1.0d0) then
            adjusted_p_recruitment = 1.0d0
        end if

        ! Probabilistically recruit an enzyme
        if (rand() < adjusted_p_recruitment) then
            ! Apply changes with probability p_change
            if (rand() < p_change) then
            nth_position = position + nth_neighbor
            if (nth_position <= size(this%histones)) then
                character(1) :: nth_histone

                nth_histone = this%histones(nth_position)

                if (current_histone == 'A' .and. nth_histone == 'U') then
                this%histones(nth_position) = 'A'
                enzyme_changes = enzyme_changes + 1
                elseif (current_histone == 'A' .and. nth_histone == 'M') then
                this%histones(nth_position) = 'U'
                enzyme_changes = enzyme_changes + 1
                elseif (current_histone == 'M' .and. nth_histone == 'U') then
                this%histones(nth_position) = 'M'
                enzyme_changes = enzyme_changes + 1
                elseif (current_histone == 'M' .and. nth_histone == 'A') then
                this%histones(nth_position) = 'U'
                enzyme_changes = enzyme_changes + 1
                end if
            end if
            end if
        end if
        end if
    end subroutine Chromatine_ChangeNextHistones

        subroutine Polymerase_Move(this)
      class(Polymerase), intent(inout) :: this

      ! Move the polymerase based on probabilities
      if (rand() < left_movement_probability) then
        this%position = max(1, this%position - 1)
      else
        this%position = min(size(this%chromatine%histones), this%position + 1)
      end if
    end subroutine Polymerase_Move

    subroutine Polymerase_ChangeHistones(this, p_change, noisy_changes)
      class(Polymerase), intent(inout) :: this
      real(8), intent(in) :: p_change
      integer, intent(inout) :: noisy_changes

      ! Change histones in the chromatine based on probabilities
      call this%chromatine%change_next_histones(this%position, p_change, noisy_changes, 1)
    end subroutine Polymerase_ChangeHistones

    subroutine Polymerase_Delete(this, existing_polymerase_positions)
      class(Polymerase), intent(inout) :: this
      integer, dimension(:), intent(inout) :: existing_polymerase_positions

      ! Delete the polymerase by updating existing positions
      existing_polymerase_positions = pack(existing_polymerase_positions /= this%position, existing_polymerase_positions)
    end subroutine Polymerase_Delete

  end subroutine Chromatine_ChangeNextHistones

end program ChromatinSimulation


