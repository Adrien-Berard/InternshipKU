program ChromatinSimulation
  implicit none
  real(8), parameter :: alpha = 3.0d0 / 4.0d0
  real(8), parameter :: change_probability = alpha
  real(8), parameter :: regeneration_probability = 0.3d0
  real(8), parameter :: adding_polymerase_probability = 0.3d0
  real(8), parameter :: noisy_transition_probability = 1.0d0 - alpha
  real(8), parameter :: vicinity_size = 5.0d0
  real(8), parameter :: F = alpha / (1.0d0 - alpha)
  real(8), parameter :: slope = 1.0d-5
  real(8), parameter :: intercept = 0.0d0
  real(8), parameter :: left_movement_probability = 0.5d0
  real(8), parameter :: right_movement_probability = 0.5d0
  integer, parameter :: chromatine_size = 60
  integer, parameter :: polymerase_count = 0
  integer, parameter :: simulation_steps = 5000
  integer, parameter :: adding_position = 25
  integer, parameter :: end_of_replication_position = chromatine_size - 25

  integer :: histone_modification_percentage
  real(8) :: recruitment_probability
  integer :: noisy_changes, enzyme_changes
  real(8) :: p_recruitment, p_change
  integer :: position, nth_neighbor

  real(8), dimension(:), allocatable :: histones
  integer, dimension(:), allocatable :: existing_polymerase_positions
  integer, dimension(:), allocatable :: polymerase_count_over_time
  integer, dimension(:), allocatable :: active_histone_count_over_time
  integer, dimension(:), allocatable :: acetylated_histone_count_over_time
  integer, dimension(:), allocatable :: methylated_histone_count_over_time
  integer, dimension(:), allocatable :: unmodified_histone_count_over_time
  integer, dimension(:), allocatable :: noisy_changes_count_over_time
  integer, dimension(:), allocatable :: enzyme_changes_count_over_time

  integer :: frame, i
  real(8) :: local_density, probability

  ! Initialize simulation parameters
  histone_modification_percentage = 50
  recruitment_probability = 1.0d0

  ! Set seed for reproducibility
  call random_seed()
  call random_seed(size = i)
  do i = 1, i
     call random_number()
  end do

  ! Initialize arrays
  allocate(histones(chromatine_size))
  allocate(existing_polymerase_positions(0:polymerase_count))
  allocate(polymerase_count_over_time(0:simulation_steps))
  allocate(active_histone_count_over_time(0:simulation_steps))
  allocate(acetylated_histone_count_over_time(0:simulation_steps))
  allocate(methylated_histone_count_over_time(0:simulation_steps))
  allocate(unmodified_histone_count_over_time(0:simulation_steps))
  allocate(noisy_changes_count_over_time(0:simulation_steps))
  allocate(enzyme_changes_count_over_time(0:simulation_steps))

  ! Initialize chromatine and polymerases
  call initialize_chromatine(histones, histone_modification_percentage)
  call initialize_polymerases(existing_polymerase_positions, polymerase_count, adding_position)

  ! Main simulation loop
  do frame = 1, simulation_steps
     call update(frame, histones, existing_polymerase_positions, noisy_changes, enzyme_changes, &
                  polymerase_count_over_time, active_histone_count_over_time, acetylated_histone_count_over_time, &
                  methylated_histone_count_over_time, unmodified_histone_count_over_time, &
                  noisy_changes_count_over_time, enzyme_changes_count_over_time)
  end do

  ! Deallocate arrays
  deallocate(histones, existing_polymerase_positions, polymerase_count_over_time, active_histone_count_over_time, &
             acetylated_histone_count_over_time, methylated_histone_count_over_time, unmodified_histone_count_over_time, &
             noisy_changes_count_over_time, enzyme_changes_count_over_time)

contains

  subroutine initialize_chromatine(histones, histone_modification_percentage)
    real(8), intent(out), dimension(:) :: histones
    integer, intent(in) :: histone_modification_percentage

    ! Initialize chromatine with histones
    histones = 'M'
    histones(1:histone_modification_percentage) = 'U'

    ! Randomly set approximately histone_modification_percentage of histones to 'U' (unmodified)
    call random_shuffle(histones)
  end subroutine initialize_chromatine

  subroutine initialize_polymerases(existing_polymerase_positions, polymerase_count, adding_position)
    integer, intent(out), dimension(:) :: existing_polymerase_positions
    integer, intent(out) :: polymerase_count
    integer, intent(in) :: adding_position
    integer :: i, new_position

    ! Initialize polymerases with non-overlapping positions
    polymerase_count = 0
    do i = 1, 10 ! Assuming a fixed number of initial polymerases (you can adjust as needed)
       new_position = adding_position + i - 1
       existing_polymerase_positions(i) = new_position
       polymerase_count = polymerase_count + 1
    end do
  end subroutine initialize_polymerases

  subroutine update(frame, histones, existing_polymerase_positions, noisy_changes, enzyme_changes, &
                    polymerase_count_over_time, active_histone_count_over_time, acetylated_histone_count_over_time, &
                    methylated_histone_count_over_time, unmodified_histone_count_over_time, &
                    noisy_changes_count_over_time, enzyme_changes_count_over_time)
    integer, intent(in) :: frame
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(inout), dimension(:) :: existing_polymerase_positions
    integer, intent(inout) :: noisy_changes, enzyme_changes
    integer, intent(inout), dimension(:) :: polymerase_count_over_time
    integer, intent(inout), dimension(:) :: active_histone_count_over_time
    integer, intent(inout), dimension(:) :: acetylated_histone_count_over_time
    integer, intent(inout), dimension(:) :: methylated_histone_count_over_time
    integer, intent(inout), dimension(:) :: unmodified_histone_count_over_time
    integer, intent(inout), dimension(:) :: noisy_changes_count_over_time
    integer, intent(inout), dimension(:) :: enzyme_changes_count_over_time
    real(8) :: local_density, probability
    integer :: position, nth_neighbor
    integer, dimension(:), allocatable :: deleted_positions, polymerase_positions
    integer :: i, new_polymerase_positions

    ! Allocate temporary arrays
    allocate(deleted_positions(0:size(existing_polymerase_positions)))
    allocate(polymerase_positions(0:size(existing_polymerase_positions)))

    ! Update polymerase positions and histones
    do i = 1, size(existing_polymerase_positions)
       call move_polymerase(existing_polymerase_positions(i), polymerase_positions, deleted_positions)
       call change_histones(histones, existing_polymerase_positions(i))
    end do

    ! Change the next histones based on the influence of first neighbors
    position = 1 + int(random_number() * (chromatine_size - 1))
    nth_neighbor = 1 + int(random_number() * (chromatine_size - 1))
    call change_next_histones(histones, position, recruitment_probability, change_probability, enzyme_changes, nth_neighbor)
    call noisy_transition(histones, position, noisy_transition_probability, noisy_changes)

    ! Regenerate histones at unmodified positions (commented out as per original Python code)
    ! if (random_number() < regeneration_probability) then
    !    call regenerate_histones(histones, deleted_positions)
    ! end if

    ! Randomly add new polymerase at the beginning of the chromatine with a certain probability
    if (random_number() < adding_polymerase_probability) then
       call add_polymerases(existing_polymerase_positions, adding_position, new_polymerase_positions)
       call initialize_new_polymerases(histones, new_polymerase_positions)
    end if

    ! Update the number of polymerases and active histones lists
    polymerase_count_over_time(frame + 1) = size(existing_polymerase_positions)
    active_histone_count_over_time(frame + 1) = count(histones == 'M' .or. histones == 'A')
    acetylated_histone_count_over_time(frame + 1) = count(histones == 'A')
    methylated_histone_count_over_time(frame + 1) = count(histones == 'M')
    unmodified_histone_count_over_time(frame + 1) = count(histones == 'U')
    noisy_changes_count_over_time(frame + 1) = noisy_changes
    enzyme_changes_count_over_time(frame + 1) = enzyme_changes

    ! Clear the previous frame after updating the data
    call clear_previous_frame()
  end subroutine update

  subroutine move_polymerase(position, polymerase_positions, deleted_positions)
    integer, intent(in) :: position
    integer, intent(inout), dimension(:) :: polymerase_positions
    integer, intent(inout), dimension(:) :: deleted_positions
    integer :: next_position, i, total_polymerases

    ! Define two possible states for movement (left and right)
    total_polymerases = size(polymerase_positions)
    next_position = position + 1

    ! Check if there is another polymerase in the next position
    do i = 1, total_polymerases
       if (next_position == polymerase_positions(i)) then
          ! Do not move if there is another polymerase 1 place after
          return
       end if
    end do

    ! Update polymerase position
    call choose_next_position(position, next_position, polymerase_positions)

    ! Bounding conditions
    if (next_position == end_of_replication_position) then
       call delete_polymerase(position, polymerase_positions, deleted_positions)
    end if
  end subroutine move_polymerase

  subroutine choose_next_position(position, next_position, polymerase_positions)
    integer, intent(in) :: position, next_position
    integer, intent(inout), dimension(:) :: polymerase_positions
    real(8), dimension(2) :: probabilities
    real(8) :: total_prob, normalized_probabilities(2)

    ! Check if next position is already occupied by another polymerase
    if (any(next_position == polymerase_positions)) return

    ! Define two possible states for movement (left and right)
    probabilities = [left_movement_probability, right_movement_probability]

    ! Normalize probabilities to sum to 1
    total_prob = sum(probabilities)
    normalized_probabilities = probabilities / total_prob

    ! Choose the next position based on normalized probabilities
    if (random_number() < normalized_probabilities(1)) then
       polymerase_positions(position) = next_position
    end if
  end subroutine choose_next_position

  subroutine delete_polymerase(position, polymerase_positions, deleted_positions)
    integer, intent(in) :: position
    integer, intent(inout), dimension(:) :: polymerase_positions
    integer, intent(inout), dimension(:) :: deleted_positions
    integer :: i

    ! Remove polymerase from polymerase_positions
    do i = 1, size(polymerase_positions)
       if (position == polymerase_positions(i)) then
          deleted_positions(i) = position
          polymerase_positions(i) = 0
          exit
       end if
    end do
  end subroutine delete_polymerase

  subroutine change_histones(histones, position)
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(in) :: position
    character(1) :: current_histone, nth_histone

    ! Simulate the histone change process by polymerase
    if (position >= 1 .and. position <= size(histones)) then
       current_histone = histones(position)

       ! Calculate the influence of vicinity on the recruitment probability
       call calculate_recruitment_influence(position, recruitment_probability)

       ! Probabilistically recruit an enzyme
       if (random_number() < recruitment_probability) then
          ! Apply changes with probability change_probability
          if (random_number() < change_probability) then
             call apply_changes(histones, position)
          end if
       end if
    end if
  end subroutine change_histones

  subroutine calculate_recruitment_influence(position, recruitment_probability)
    integer, intent(in) :: position
    real(8), intent(inout) :: recruitment_probability
    real(8) :: adjusted_p_recruitment

    ! Calculate the influence of vicinity on the recruitment probability
    adjusted_p_recruitment = recruitment_probability
    ! adjusted_p_recruitment = recruitment_probability / (1.0d-1 * nth_neighbor) ! Add 1 to avoid division by zero
    if (adjusted_p_recruitment > 1.0d0) adjusted_p_recruitment = 1.0d0

    ! TODO: Implement the rest of the calculation if needed
  end subroutine calculate_recruitment_influence

  subroutine apply_changes(histones, position)
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(in) :: position
    character(1) :: current_histone, nth_histone

    ! Apply changes based on the current and neighboring histones
    if (position >= 1 .and. position < size(histones) - 1) then
       current_histone = histones(position)

       ! TODO: Implement the rest of the changes based on the original Python code
    end if
  end subroutine apply_changes

  subroutine change_next_histones(histones, position, p_recruitment, p_change, enzyme_changes, nth_neighbor)
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(in) :: position, nth_neighbor
    real(8), intent(in) :: p_recruitment, p_change
    integer, intent(inout) :: enzyme_changes
    character(1) :: current_histone, nth_histone
    real(8) :: adjusted_p_recruitment

    ! Simulate the influence of neighbors on the next histones
    if (position >= 1 .and. position < size(histones) - 1) then
       current_histone = histones(position)

       ! Calculate the influence of vicinity on the recruitment probability
       call calculate_recruitment_influence(position, recruitment_probability)

       ! Probabilistically recruit an enzyme
       if (random_number() < adjusted_p_recruitment) then
          ! Apply changes with probability p_change
          if (random_number() < p_change) then
             nth_position = position + nth_neighbor
             if (nth_position <= size(histones)) then
                nth_histone = histones(nth_position)

                ! TODO: Implement the rest of the changes based on the original Python code
             end if
          end if
       end if
    end if
  end subroutine change_next_histones

  subroutine noisy_transition(histones, position, noisy_transition_probability, noisy_changes)
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(in) :: position
    real(8), intent(in) :: noisy_transition_probability
    integer, intent(inout) :: noisy_changes

    ! Simulate noisy transitions
    if (random_number() < noisy_transition_probability / 3.0d0) then
       if (histones(position) == 'A') then
          histones(position) = 'U'
          noisy_changes = noisy_changes + 1
       elseif (histones(position) == 'M') then
          histones(position) = 'U'
          noisy_changes = noisy_changes + 1
       elseif (histones(position) == 'U') then
          if (random_number() < 0.5d0) then
             histones(position) = 'A'
             noisy_changes = noisy_changes + 1
          else
             histones(position) = 'M'
             noisy_changes = noisy_changes + 1
          end if
       end if
    end if
  end subroutine noisy_transition

  subroutine regenerate_histones(histones, deleted_positions)
    real(8), intent(inout), dimension(:) :: histones
    integer, intent(in), dimension(:) :: deleted_positions
    integer :: i

    ! Regenerate histones at unmodified positions
    do i = 1, size(deleted_positions)
       if (histones(deleted_positions(i)) == 'U') then
          ! TODO: Implement the regeneration process based on the original Python code
       end if
    end do
  end subroutine regenerate_histones

  subroutine add_polymerases(existing_polymerase_positions, adding_position, new_polymerase_positions)
    integer, intent(inout), dimension(:) :: existing_polymerase_positions
    integer, intent(in) :: adding_position
    integer, intent(out), dimension(:) :: new_polymerase_positions
    integer :: i, new_position

    ! Add a specified number of new polymerases at non-overlapping positions
    do i = 1, 1 ! Assuming a fixed number of new polymerases to be added (you can adjust as needed)
       new_position = adding_position
       do while (any(new_position == existing_polymerase_positions))
          new_position = new_position + 1
       end do
       existing_polymerase_positions(size(existing_polymerase_positions) + 1) = new_position
       new_polymerase_positions(i) = new_position
    end do
  end subroutine add_polymerases

  subroutine initialize_new_polymerases(histones, new_polymerase_positions)
    real(8), intent(in) :: histones
    integer, intent(in), dimension(:) :: new_polymerase_positions
    integer :: i

    ! Initialize new polymerases with specified positions
    do i = 1, size(new_polymerase_positions)
       call initialize_polymerase(histones, new_polymerase_positions(i))
    end do
  end subroutine initialize_new_polymerases

  subroutine initialize_polymerase(histones, position)
    real(8), intent(in) :: histones
    integer, intent(in) :: position

    ! TODO: Implement the initialization of a new polymerase based on the original Python code
  end subroutine initialize_polymerase

  subroutine clear_previous_frame()
    ! Clear the previous frame in a way suitable for your Fortran environment
    ! This may involve clearing the console or any other output method used in your environment
    ! TODO: Implement based on your Fortran environment
  end subroutine clear_previous_frame

end program ChromatinSimulation
