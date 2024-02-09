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
  real(8), parameter :: this

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

  ! Initialize chromatine and polymerases with a specified temperature
  call chromatine%initialize(chromatine_size)
  call initialize_polymerases(existing_polymerase_positions, polymerase_count, adding_position)

  ! Create an animated plot
  call AnimatedPlot()

  contains

  ! Subroutine to initialize chromatine
  subroutine Chromatine_Initialize(this, histones_count)
    class(Chromatine), intent(inout) :: this
    integer, intent(in) :: histones_count

    ! Initialize chromatine with histones
    allocate(this%histones(histones_count))
    this%histones = 'M'

    ! Randomly set approximately histone_modification_percentage of histones to 'U' (unmodified)
    this%histones(1:int(histone_modification_percentage * histones_count)) = 'U'
  end subroutine Chromatine_Initialize

  ! Subroutine to simulate noisy transition in chromatine
  subroutine Chromatine_NoisyTransition(this, position, noisy_transition_probability, noisy_changes)
    class(Chromatine), intent(inout) :: this
    integer, intent(in) :: position
    real(8), intent(in) :: noisy_transition_probability
    integer, intent(inout) :: noisy_changes

    if (rand() < noisy_transition_probability / 3.0d0) then
      if (this%histones(position) == 'A' .or. this%histones(position) == 'M') then
        this%histones(position) = 'U'
        noisy_changes = noisy_changes + 1
      elseif (this%histones(position) == 'U') then
        if (rand() < 0.5d0) then
          this%histones(position) = 'A'
        else
          this%histones(position) = 'M'
        end if
        noisy_changes = noisy_changes + 1
      end if
    end if
  end subroutine Chromatine_NoisyTransition

  ! Subroutine to add polymerases to chromatine
  subroutine Chromatine_AddPolymerases(this, count, existing_polymerase_positions, adding_position)
    class(Chromatine), intent(inout) :: this
    integer, intent(in) :: count
    integer, dimension(:), intent(inout) :: existing_polymerase_positions
    integer, intent(in) :: adding_position

    integer :: i, new_position

    ! Add a specified number of new polymerases at non-overlapping positions
    do i = 1, count
      new_position = adding_position
      do while (new_position `in` existing_polymerase_positions)
        new_position = new_position + 1
      end do
      existing_polymerase_positions = [existing_polymerase_positions, new_position]
    end do
  end subroutine Chromatine_AddPolymerases

  ! Function to calculate the probability of adding a new polymerase
  real(8) function Chromatine_AddingPolyProba(this, adding_position)
    class(Chromatine), intent(in) :: this
    integer, intent(in) :: adding_position

    integer :: start_index, end_index, local_density

    ! Linear function of the local density of histones
    ! Let's calculate local density as the count of active histones in the vicinity of the polymerase
    start_index = max(1, adding_position - vicinity_size)
    end_index = min(size(this%histones), adding_position + vicinity_size + 1)

    local_density = count(this%histones(start_index:end_index) == 'M' .or. this%histones(start_index:end_index) == 'A')

    ! Calculate the probability of adding a new polymerase
    Chromatine_AddingPolyProba = slope * real(local_density) + intercept
  end function Chromatine_AddingPolyProba

  ! Subroutine to simulate the influence of neighbors on the next histones
  subroutine Chromatine_ChangeNextHistones(this, position, change_probability, noisy_changes)
    class(Chromatine), intent(inout) :: this
    integer, intent(in) :: position
    real(8), intent(in) :: change_probability
    integer, intent(inout) :: noisy_changes

    integer :: left_position, right_position

    ! Calculate left and right neighbor positions
    left_position = max(1, position - 1)
    right_position = min(size(this%histones), position + 1)

    ! Check if either neighbor histone is 'U', and if yes, change the histone at the current position to 'U'
    if (this%histones(left_position) == 'U' .or. this%histones(right_position) == 'U') then
      this%histones(position) = 'U'
      noisy_changes = noisy_changes + 1
    else
      ! Check if a histone at the current position should change based on a probability
      if (rand() < change_probability) then
        this%histones(position) = 'U'
        noisy_changes = noisy_changes + 1
      end if
    end if
  end subroutine Chromatine_ChangeNextHistones

  ! Subroutine to initialize polymerases
  subroutine initialize_polymerases(existing_polymerase_positions, count, adding_position)
    integer, dimension(:), allocatable :: existing_polymerase_positions
    integer, intent(in) :: count
    integer, intent(in) :: adding_position

    ! Initialize existing polymerase positions
    allocate(existing_polymerase_positions(0:count))
    existing_polymerase_positions = adding_position

    ! Initialize the chromatine with polymerases
    call chromatine%add_polymerases(count, existing_polymerase_positions, adding_position)
  end subroutine initialize_polymerases

  ! Subroutine to move polymerases
  subroutine Polymerase_Move(this)
    class(Polymerase), intent(inout) :: this

    ! Move polymerase left or right based on probabilities
    if (rand() < left_movement_probability) then
      this%position = max(1, this%position - 1)
    else
      this%position = min(size(this%chromatine%histones), this%position + 1)
    end if
  end subroutine Polymerase_Move

  ! Subroutine to change histones at the polymerase position
  subroutine Polymerase_ChangeHistones(this, change_probability, noisy_changes)
    class(Polymerase), intent(inout) :: this
    real(8), intent(in) :: change_probability
    integer, intent(inout) :: noisy_changes

    ! Change the histones at the polymerase position based on a probability
    call this%chromatine%change_next_histones(this%position, change_probability, noisy_changes)
  end subroutine Polymerase_ChangeHistones

  ! Subroutine to delete polymerases
  subroutine Polymerase_Delete(this, existing_polymerase_positions)
    class(Polymerase), intent(inout) :: this
    integer, dimension(:), intent(inout) :: existing_polymerase_positions

    ! Remove the current polymerase position from the list of existing polymerase positions
    existing_polymerase_positions = existing_polymerase_positions(/1:this%position-1, this%position+1:/)
  end subroutine Polymerase_Delete

  ! Subroutine to plot the chromatine state and polymerase positions over time
  subroutine AnimatedPlot()
    integer :: frame, noisy_changes, enzyme_changes, i

    ! Initialize counts
    noisy_changes = 0
    enzyme_changes = 0

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

      ! Plot the chromatine state and polymerase positions
      call PlotChromatineAndPolymerases(chromatine, existing_polymerase_positions, frame)
    end do
  end subroutine AnimatedPlot

  ! Subroutine to plot the chromatine state and polymerase positions
  subroutine PlotChromatineAndPolymerases(chromatine, existing_polymerase_positions, frame)
    class(Chromatine), intent(in) :: chromatine
    integer, dimension(:), intent(in) :: existing_polymerase_positions
    integer, intent(in) :: frame

    ! Placeholder for plotting logic
    ! Implement your plotting code here based on the chromatine and polymerase positions
    ! You can use the frame variable to track the simulation steps
  end subroutine PlotChromatineAndPolymerases

  ! Function to generate a random number between 0 and 1
  real(8) function rand()
    call random_number(rand)
  end function rand

  ! Function to seed the random number generator
  subroutine srand(seed)
    integer, intent(in) :: seed
    call random_seed()
    call random_seed(size = i)
    call random_seed(put = seed, iostat = i)
  end subroutine srand

end program ChromatinSimulation
