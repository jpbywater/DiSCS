import random

def generate_sequence(number_of_states=2,
                      number_of_actions_in_each_state=3,
                      p_between_state=0.1,
                      sequence_length=20,
                      p_common_action = 0):

    state_names = ["A", "B", "C", "D", "E", "F", "G"][:number_of_states]
    state_id = 0
    action_id = 0
    output_state_names = []
    output_action_names = []

    for i in range(sequence_length):
        # save current state and action to output lists
        output_state_names.append(state_names[state_id])
        if random.random() < p_common_action:
            output_action_names.append('common')
        else:
            output_action_names.append(state_names[state_id].lower() + str(action_id))
        # transition state
        if random.random() < p_between_state:
            # transition to another state randomly (equally likely excluding current).
            state_ids = list(range(number_of_states))
            state_ids.remove(state_id)
            state_id = random.choice(state_ids)
        # whether transition state or not, randomly select new action id (equally likely including current)
        action_ids=list(range(number_of_actions_in_each_state))
        action_id=random.choice(action_ids)
    return output_state_names, output_action_names

